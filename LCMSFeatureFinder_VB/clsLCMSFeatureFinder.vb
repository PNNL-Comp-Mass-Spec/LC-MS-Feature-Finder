Option Strict On

Imports System.IO
Imports System.Reflection
Imports System.Runtime.InteropServices
Imports System.Threading
Imports System.Windows.Forms
Imports UMCCreation
' This class reads a text file with mass and intensity data for MS spectra
' and determines the LC-MS features present using UMCCreation.dll
'
' -------------------------------------------------------------------------------
' Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
' Program started January 30, 2007
' Copyright 2007, Battelle Memorial Institute.  All Rights Reserved.

' E-mail: matthew.monroe@pnl.gov or matt@alchemistmatt.com
' Website: http://omics.pnl.gov/ or http://www.sysbio.org/resources/staff/ or http://panomics.pnnl.gov/
' -------------------------------------------------------------------------------
'
' Licensed under the Apache License, Version 2.0; you may not use this file except
' in compliance with the License.  You may obtain a copy of the License at
' http://www.apache.org/licenses/LICENSE-2.0
'
' Notice: This computer software was prepared by Battelle Memorial Institute,
' hereinafter the Contractor, under Contract No. DE-AC05-76RL0 1830 with the
' Department of Energy (DOE).  All rights in the computer software are reserved
' by DOE on behalf of the United States Government and the Contractor as
' provided in the Contract.  NEITHER THE GOVERNMENT NOR THE CONTRACTOR MAKES ANY
' WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS
' SOFTWARE.  This notice including this sentence must appear on any copies of
' this computer software.

Public Class clsLCMSFeatureFinder
    Public Event ProgressReset()
    Public Event ProgressChanged(taskDescription As String, percentComplete As Single)     ' PercentComplete ranges from 0 to 100, but can contain decimal percentage values
    Public Event ProgressComplete()

    Public Sub New()
        mFileDate = PROGRAM_DATE
        InitializeLocalVariables()
    End Sub

#Region "Constants and Enums"
    Protected Const COLUMN_DELIMITER As Char = ControlChars.Tab
    Protected Const FILE_SUFFIX_LCMS_FEATURES As String = "_Features"
    Protected Const FILE_SUFFIX_PEAK_TO_FEATURE_MAP As String = "_PeakToFeatureMap"

    Public Enum eLCMSFeatureFinderErrorCodes
        NoError = 0
        InvalidInputFilePath = 1
        InvalidOutputFolderPath = 2
        ParameterFileNotFound = 4
        FilePathError = 8

        ParameterFileReadError = 16
        UnknownFileExtension = 32
        InputFileReadError = 64
        OutputFileWriteError = 128
        LCMSProcessingError = 256

        UnspecifiedError = -1
    End Enum

    Protected Const INPUT_FILE_COLUMN_NAME_COUNT As Integer = 15
    Protected Enum eInputFileColumnNames
        Scan = 0
        Charge = 1
        Abundance = 2
        MZOfMostAbu = 3
        Fit = 4
        AverageMass = 5
        MonoisotopicMass = 6
        MassOfMostAbundant = 7
        FWHM = 8
        SignalToNoise = 9
        MonoMassAbundance = 10
        MonoMassPlus2DaAbundance = 11
        Index = 12
        NET = 13
        IMSDriftTime = 14
    End Enum

    ' Note: These should all be lowercase string values
    Protected Const ISOS_COLUMN_SCAN_NUM As String = "scan_num"
    Protected Const ISOS_COLUMN_FRAME_NUM As String = "frame_num"
    Protected Const ISOS_COLUMN_IMS_SCAN_NUM As String = "ims_scan_num"         ' LC-MS Feature Finder ignores this column
    Protected Const ISOS_COLUMN_CHARGE As String = "charge"
    Protected Const ISOS_COLUMN_ABUNDANCE As String = "abundance"
    Protected Const ISOS_COLUMN_MZ As String = "mz"
    Protected Const ISOS_COLUMN_FIT As String = "fit"
    Protected Const ISOS_COLUMN_AVERAGE_MW As String = "average_mw"
    Protected Const ISOS_COLUMN_MONOISOTOPIC_MW As String = "monoisotopic_mw"
    Protected Const ISOS_COLUMN_MOSTABUNDANT_MW As String = "mostabundant_mw"
    Protected Const ISOS_COLUMN_FWHM As String = "fwhm"
    Protected Const ISOS_COLUMN_SIGNAL_NOISE As String = "signal_noise"
    Protected Const ISOS_COLUMN_MONO_ABUNDANCE As String = "mono_abundance"
    Protected Const ISOS_COLUMN_MONO_PLUS2_ABUNDANCE As String = "mono_plus2_abundance"
    Protected Const ISOS_COLUMN_INDEX As String = "index"                       ' Index of the data point in the source application
    Protected Const ISOS_COLUMN_NET As String = "net"                           ' Normalized elution time; not used by LC-MS Feature Finder
    Protected Const ISOS_COLUMN_IMS_DRIFT_TIME As String = "ims_drift_time"     ' IMS Drift Time column in Decon2LS v1
    Protected Const ISOS_COLUMN_IMS_DRIFT_TIME_ALT As String = "drift_time"     ' IMS Drift Time column in Decon2LS v2
    Protected Const ISOS_COLUMN_ORIG_INTENSITY As String = "orig_intensity"     ' LC-MS Feature Finder ignores this column
    Protected Const ISOS_COLUMN_TIA_ORIG_INTENSITY As String = "tia_orig_intensity"     ' LC-MS Feature Finder ignores this column

    Protected Const INI_SECTION_UMC_CREATION_OPTIONS As String = "UMCCreationOptions"
    Protected Enum eIniFileSectionConstants
        UnknownSection = 0
        UMCCreationOptions = 1
    End Enum
#End Region

#Region "Structures"
    ' Set any of the weights to 0 to disable their use
    Public Structure udtFeatureFindingOptionsType
        Public MonoMassWeight As Single
        Public MonoMassConstraint As Single
        Public MonoMassConstraintIsPPM As Boolean

        Public AvgMassWeight As Single
        Public AvgMassConstraint As Single
        Public AvgMassConstraintIsPPM As Boolean

        Public LogAbundanceWeight As Single
        Public ScanWeight As Single         ' Used if .UseGenericNET is False
        Public NETWeight As Single          ' Used if .UseGenericNET is True
        Public FitWeight As Single
        Public IMSDriftTimeWeight As Single

        Public MaxDistance As Single
        Public UseGenericNET As Boolean

        Public MinScan As Integer
        Public MaxScan As Integer
        Public MinFeatureLengthPoints As Integer
        Public RequireMatchingChargeState As Boolean   ' When true, then all data points in the UMC will have the same charge state
    End Structure

#End Region

#Region "Classwide Variables"
    Protected mFileDate As String
    Protected mErrorCode As eLCMSFeatureFinderErrorCodes

    Protected mProgressStepDescription As String
    Protected mProgressPercentComplete As Single        ' Ranges from 0 to 100, but can contain decimal percentage values
    Protected mStatusMessage As String
    Protected mAbortProcessing As Boolean

    Protected mIsotopePeaksCount As Integer
    Protected mIsotopePeaks() As clsIsotopePeak
    Protected mIsotopePeaksMapToIndexInDataSource() As Integer      ' Parallel to mIsotopePeaks

    Protected mUMCCreator As clsUMCCreator

    ' These values are useful when calling this program with just a narrow range of data points.  We need to know
    '  the min and max scan range so that we can convert from scan to NET, and if we only have a narrow range of points
    '  then the conversion to NET will be incorrect
    ' If both of these values are 0, then the min and max scan values are determined from the data
    Public mFeatureFindingOptions As udtFeatureFindingOptionsType

#End Region

#Region "Properties"
    Public Property AbortProcessing() As Boolean
        Get
            Return mAbortProcessing
        End Get
        Set(Value As Boolean)
            mAbortProcessing = Value
        End Set
    End Property

    Public ReadOnly Property ErrorCode() As eLCMSFeatureFinderErrorCodes
        Get
            Return mErrorCode
        End Get
    End Property

    Public Overridable ReadOnly Property ProgressStepDescription() As String
        Get
            Return mProgressStepDescription
        End Get
    End Property

    ' ProgressPercentComplete ranges from 0 to 100, but can contain decimal percentage values
    Public ReadOnly Property ProgressPercentComplete() As Single
        Get
            Return CType(Math.Round(mProgressPercentComplete, 2), Single)
        End Get
    End Property

    Public ReadOnly Property StatusMessage() As String
        Get
            Return mStatusMessage
        End Get
    End Property
#End Region

    Public Sub AddIsotopePeak(intScanNumber As Integer,
                                    intIndexInDataSource As Integer,
                                    dblAbundance As Double,
                                    dblMonoMassPlus2DaAbundance As Double,
                                    dblMonoisotopicMass As Double,
                                    dblAverageMass As Double,
                                    dblMassOfMostAbundant As Double,
                                    dblMZOfMostAbu As Double,
                                    intCharge As Short,
                                    sngFit As Single,
                                    sngIMSDriftTime As Single)

        If mIsotopePeaks Is Nothing Then ClearIsotopePeaks()

        If mIsotopePeaksCount >= mIsotopePeaks.Length Then
            ReDim Preserve mIsotopePeaks(mIsotopePeaksCount * 2 - 1)
            ReDim Preserve mIsotopePeaksMapToIndexInDataSource(mIsotopePeaksCount * 2 - 1)
        End If

        Dim originalIndex = mIsotopePeaksCount

        mIsotopePeaks(originalIndex) = New clsIsotopePeak(originalIndex)
        With mIsotopePeaks(mIsotopePeaksCount)
            .mint_lc_scan = intScanNumber

            ' Note: the originalIndex value must be 0 for the first isotope peak, and then increment up from there
            ' The following array is used to track the intIndexInDataSource value for each entry in mIsotopePeaks
            mIsotopePeaksMapToIndexInDataSource(originalIndex) = intIndexInDataSource

            .mdbl_abundance = dblAbundance
            ' .mdbl_mono_abundance = GetColumnValueDbl(strSplitLine, intColumnMapping(eInputFileColumnNames.MonoMassAbundance), 0)
            .mdbl_i2_abundance = dblMonoMassPlus2DaAbundance

            .mdbl_mono_mass = dblMonoisotopicMass
            .mdbl_average_mass = dblAverageMass
            .mdbl_max_abundance_mass = dblMassOfMostAbundant
            .mdbl_mz = dblMZOfMostAbu
            .mshort_charge = intCharge

            .mflt_fit = sngFit
            .mflt_ims_drift_time = sngIMSDriftTime
        End With
        mIsotopePeaksCount += 1

    End Sub

    Protected Sub AutoLoadOptions(strInputFilePath As String)
        ' ToDo: Look for a .Ini file that matches strInputFilePath
        ' Open it and parse the settings present, updatingmFeatureFindingOptions

        Dim strIniFilePath As String = String.Empty
        Dim strInputFolderPath As String

        Dim objIniFileReader As clsIniFileReader

        Try
            If strInputFilePath Is Nothing OrElse strInputFilePath.Length = 0 Then
                Exit Sub
            End If

            strInputFolderPath = GetParentFolderPathForFile(strInputFilePath)
            strIniFilePath = Path.Combine(strInputFolderPath, Path.GetFileNameWithoutExtension(strInputFilePath) & ".ini")

            If Not File.Exists(strIniFilePath) Then
                Console.WriteLine("Settings file not found (" & strIniFilePath & "); using defaults")
                Return
            End If

            Console.WriteLine("Reading settings from " & strIniFilePath)

            ' Read the file
            objIniFileReader = New clsIniFileReader(strIniFilePath)

            If Not objIniFileReader.DataLoaded Then
                Console.WriteLine("Settings were not loaded from the settings file (" & strIniFilePath & "); likely an invalid format")
                Return
            End If

            With mFeatureFindingOptions
                .MonoMassWeight = objIniFileReader.GetSetting(INI_SECTION_UMC_CREATION_OPTIONS, "MonoMassWeight", .MonoMassWeight)
                .MonoMassConstraint = objIniFileReader.GetSetting(INI_SECTION_UMC_CREATION_OPTIONS, "MonoMassConstraint", .MonoMassConstraint)
                .MonoMassConstraintIsPPM = objIniFileReader.GetSetting(INI_SECTION_UMC_CREATION_OPTIONS, "MonoMassConstraintIsPPM", .MonoMassConstraintIsPPM)

                .AvgMassWeight = objIniFileReader.GetSetting(INI_SECTION_UMC_CREATION_OPTIONS, "AvgMassWeight", .AvgMassWeight)
                .AvgMassConstraint = objIniFileReader.GetSetting(INI_SECTION_UMC_CREATION_OPTIONS, "AvgMassConstraint", .AvgMassConstraint)
                .AvgMassConstraintIsPPM = objIniFileReader.GetSetting(INI_SECTION_UMC_CREATION_OPTIONS, "AvgMassConstraintIsPPM", .AvgMassConstraintIsPPM)

                .LogAbundanceWeight = objIniFileReader.GetSetting(INI_SECTION_UMC_CREATION_OPTIONS, "LogAbundanceWeight", .LogAbundanceWeight)
                .ScanWeight = objIniFileReader.GetSetting(INI_SECTION_UMC_CREATION_OPTIONS, "ScanWeight", .ScanWeight)
                .NETWeight = objIniFileReader.GetSetting(INI_SECTION_UMC_CREATION_OPTIONS, "NETWeight", .NETWeight)
                .FitWeight = objIniFileReader.GetSetting(INI_SECTION_UMC_CREATION_OPTIONS, "FitWeight", .FitWeight)
                .IMSDriftTimeWeight = objIniFileReader.GetSetting(INI_SECTION_UMC_CREATION_OPTIONS, "IMSDriftTimeWeight", .IMSDriftTimeWeight)

                .MaxDistance = objIniFileReader.GetSetting(INI_SECTION_UMC_CREATION_OPTIONS, "MaxDistance", .MaxDistance)
                .UseGenericNET = objIniFileReader.GetSetting(INI_SECTION_UMC_CREATION_OPTIONS, "UseGenericNET", .UseGenericNET)

                .MinScan = objIniFileReader.GetSetting(INI_SECTION_UMC_CREATION_OPTIONS, "MinScan", .MinScan)
                .MaxScan = objIniFileReader.GetSetting(INI_SECTION_UMC_CREATION_OPTIONS, "MaxScan", .MaxScan)

                .MinScan = objIniFileReader.GetSetting(INI_SECTION_UMC_CREATION_OPTIONS, "LCMinScan", .MinScan)
                .MaxScan = objIniFileReader.GetSetting(INI_SECTION_UMC_CREATION_OPTIONS, "LCMaxScan", .MaxScan)

                .MinFeatureLengthPoints = objIniFileReader.GetSetting(INI_SECTION_UMC_CREATION_OPTIONS, "MinFeatureLengthPoints", .MinFeatureLengthPoints)

                .RequireMatchingChargeState = objIniFileReader.GetSetting(INI_SECTION_UMC_CREATION_OPTIONS, "RequireMatchingChargeState", .RequireMatchingChargeState)
                .RequireMatchingChargeState = objIniFileReader.GetSetting(INI_SECTION_UMC_CREATION_OPTIONS, "UseCharge", .RequireMatchingChargeState)
            End With


        Catch ex As Exception
            LogErrors("AutoLoadOptions", "Error reading the .Ini file(" & Path.GetFileName(strIniFilePath) & ")", ex, False, False, True, eLCMSFeatureFinderErrorCodes.InputFileReadError)
        End Try

    End Sub

    Public Sub ClearIsotopePeaks(Optional ByVal intInitialPeakReserveCount As Integer = 1000)
        mIsotopePeaksCount = 0
        ReDim mIsotopePeaks(intInitialPeakReserveCount - 1)
        ReDim mIsotopePeaksMapToIndexInDataSource(intInitialPeakReserveCount - 1)
    End Sub

    Protected Function FindLCMSFeaturesFromFile(strInputFilePath As String, strFeaturesOutputFilePath As String, strFeatureToPeakMapFilePath As String) As Boolean
        ' Returns True if Success, False if failure
        ' Note: This function assumes strInputFilePath exists

        Const MIN_COLUMN_COUNT As Integer = 3

        Dim cColumnDelimiter As Char

        Dim ioFile As FileInfo
        Dim srOutFileFeatures As StreamWriter
        Dim srOutFileFeatureToIsotopePeakMap As StreamWriter

        Dim strMessage As String
        Dim strLineIn As String
        Dim strLineout As String
        Dim strSplitLine() As String

        Dim intDataLinesRead As Integer
        Dim intIndex As Integer

        Dim intScanNumber As Integer
        Dim intIndexInDataSource As Integer
        Dim dblAbundance As Double
        Dim dblMonoMassPlus2DaAbundance As Double
        Dim dblMonoisotopicMass As Double
        Dim dblAverageMass As Double
        Dim dblMassOfMostAbundant As Double
        Dim dblMZOfMostAbu As Double
        Dim intCharge As Short
        Dim sngFit As Single
        Dim sngIMSDriftTime As Single

        Dim blnColumnMappingDefined As Boolean
        Dim intColumnMapping() As Integer = Nothing

        Dim blnSuccess As Boolean

        ' Assume success for now
        blnSuccess = True

        Try
            strMessage = "Reading data file " & Path.GetFileName(strInputFilePath)
            Console.WriteLine(strMessage & " ")
            UpdateProgress(strMessage, 0)

            ' Obtain the full path to the file
            ioFile = New FileInfo(strInputFilePath)

            ' Open a handle to the data file
            Using srInFile = New StreamReader(New FileStream(ioFile.FullName, FileMode.Open, FileAccess.Read, FileShare.Read))

                ' Initialize the column mappings to -1
                blnColumnMappingDefined = False
                InitializeColumnMappings(intColumnMapping)

                ' Initialize define the column delimiter as a tab
                cColumnDelimiter = ControlChars.Tab
                If Path.GetExtension(strInputFilePath).ToLower = ".csv" Then
                    ' The file is a CSV file; use a comma as the delimiter
                    cColumnDelimiter = ","c
                End If

                ' Initialize the array that holds the isotope information
                ClearIsotopePeaks()

                intDataLinesRead = 0
                Do While Not srInFile.EndOfStream
                    strLineIn = srInFile.ReadLine()

                    If Not strLineIn Is Nothing AndAlso strLineIn.Length > 0 Then
                        strSplitLine = strLineIn.Split(cColumnDelimiter)

                        If Not strSplitLine Is Nothing AndAlso strSplitLine.Length >= MIN_COLUMN_COUNT Then
                            ' If this is the first row, then see if it contains column headers
                            ' If it does, then read them to update the column mapping

                            If Not blnColumnMappingDefined Then
                                ' If the first column is not a number, then define the mappings based on the column names
                                ' If it is a number, then use the default column mappings
                                If IsNumber(strSplitLine(0)) Then
                                    DefineDefaultColumnMappings(intColumnMapping)
                                Else
                                    ParseInputFileColumnNames(strSplitLine, intColumnMapping)
                                End If
                                blnColumnMappingDefined = True
                            End If

                            ' Parse strSplitLine() if the first column is a number
                            If IsNumber(strSplitLine(0)) Then
                                intScanNumber = GetColumnValueInt(strSplitLine, intColumnMapping(eInputFileColumnNames.Scan), -1)
                                If intScanNumber >= 0 Then
                                    If intColumnMapping(eInputFileColumnNames.Index) >= 0 Then
                                        intIndexInDataSource = GetColumnValueInt(strSplitLine, intColumnMapping(eInputFileColumnNames.Index), -1)
                                    Else
                                        intIndexInDataSource = intDataLinesRead
                                    End If

                                    dblAbundance = GetColumnValueDbl(strSplitLine, intColumnMapping(eInputFileColumnNames.Abundance), 0)
                                    ' .mdbl_mono_abundance = GetColumnValueDbl(strSplitLine, intColumnMapping(eInputFileColumnNames.MonoMassAbundance), 0)
                                    dblMonoMassPlus2DaAbundance = GetColumnValueDbl(strSplitLine, intColumnMapping(eInputFileColumnNames.MonoMassPlus2DaAbundance), 0)
                                    dblMonoisotopicMass = GetColumnValueDbl(strSplitLine, intColumnMapping(eInputFileColumnNames.MonoisotopicMass), 0)
                                    dblAverageMass = GetColumnValueDbl(strSplitLine, intColumnMapping(eInputFileColumnNames.AverageMass), 0)
                                    dblMassOfMostAbundant = GetColumnValueDbl(strSplitLine, intColumnMapping(eInputFileColumnNames.MassOfMostAbundant), 0)
                                    dblMZOfMostAbu = GetColumnValueSng(strSplitLine, intColumnMapping(eInputFileColumnNames.MZOfMostAbu), 0)
                                    intCharge = CShort(GetColumnValueInt(strSplitLine, intColumnMapping(eInputFileColumnNames.Charge), 1))
                                    sngFit = GetColumnValueSng(strSplitLine, intColumnMapping(eInputFileColumnNames.Fit), 0)
                                    sngIMSDriftTime = GetColumnValueSng(strSplitLine, intColumnMapping(eInputFileColumnNames.IMSDriftTime), 0)

                                    AddIsotopePeak(intScanNumber, intIndexInDataSource, dblAbundance, dblMonoMassPlus2DaAbundance,
                                                   dblMonoisotopicMass, dblAverageMass, dblMassOfMostAbundant, dblMZOfMostAbu,
                                                   intCharge, sngFit, sngIMSDriftTime)

                                    intDataLinesRead += 1
                                End If
                            End If
                        End If
                    End If
                Loop
            End Using

            ' Verify that one or more isotope peaks were loaded
            If mIsotopePeaksCount = 0 Then
                LogErrors("FindLCMSFeatures", "Did not find any valid data in the input file (" & Path.GetFileName(strInputFilePath) & ")", Nothing, True, False, True, eLCMSFeatureFinderErrorCodes.InputFileReadError)
            End If

        Catch ex As Exception
            LogErrors("FindLCMSFeatures", "Error reading the input data file (" & Path.GetFileName(strInputFilePath) & ")", ex, False, False, True, eLCMSFeatureFinderErrorCodes.InputFileReadError)
            Return False
        End Try

        Try
            ' Initialize the output files
            ioFile = New FileInfo(strFeaturesOutputFilePath)
            srOutFileFeatures = New StreamWriter(New FileStream(ioFile.FullName, FileMode.Create, FileAccess.Write, FileShare.Read))

            strLineout = "Feature_Index" & COLUMN_DELIMITER &
                         "Scan" & COLUMN_DELIMITER &
                         "Scan_Start" & COLUMN_DELIMITER &
                         "Scan_End" & COLUMN_DELIMITER &
                         "Scan_Aligned" & COLUMN_DELIMITER &
                         "NET" & COLUMN_DELIMITER &
                         "Monoisotopic_Mass" & COLUMN_DELIMITER &
                         "Monoisotopic_Mass_Calibrated" & COLUMN_DELIMITER &
                         "Abundance" & COLUMN_DELIMITER &
                         "Class_Rep_MZ" & COLUMN_DELIMITER &
                         "Class_Rep_Charge"

            srOutFileFeatures.WriteLine(strLineout)

            ioFile = New FileInfo(strFeatureToPeakMapFilePath)
            srOutFileFeatureToIsotopePeakMap = New StreamWriter(New FileStream(ioFile.FullName, FileMode.Create, FileAccess.Write, FileShare.Read))

            strLineout = "Feature_Index" & COLUMN_DELIMITER &
                         "Peak_Index"

            srOutFileFeatureToIsotopePeakMap.WriteLine(strLineout)

        Catch ex As Exception
            LogErrors("FindLCMSFeatures", "Error initializing the output file (" & Path.GetFileName(strInputFilePath) & ")", ex, False, False, True, eLCMSFeatureFinderErrorCodes.OutputFileWriteError)
            Return False
        End Try

        Dim objFeatures() As clsUMC = Nothing
        Dim intIsotopePeaksIndex() As Integer = Nothing
        Dim intUMCIndex() As Integer = Nothing

        ' Perform the work of finding the LC-MS features
        FindLCMSFeaturesFromMemory(objFeatures, intIsotopePeaksIndex, intUMCIndex)

        Try
            ' Write the Features to the output file
            For intIndex = 0 To objFeatures.Length - 1
                With objFeatures(intIndex)
                    strLineout = .mint_umc_index.ToString & COLUMN_DELIMITER &
                                 .mint_scan.ToString & COLUMN_DELIMITER &
                                 .mint_start_scan.ToString & COLUMN_DELIMITER &
                                 .mint_end_scan.ToString & COLUMN_DELIMITER &
                                 .mint_scan_aligned.ToString & COLUMN_DELIMITER &
                                 .mdbl_net.ToString("0.0000") & COLUMN_DELIMITER &
                                 ValueToString(.mdbl_mono_mass, 10) & COLUMN_DELIMITER &
                                 ValueToString(.mdbl_mono_mass_calibrated, 10) & COLUMN_DELIMITER &
                                 ValueToString(.mdbl_abundance, 8) & COLUMN_DELIMITER &
                                 ValueToString(.mdbl_class_rep_mz, 10) & COLUMN_DELIMITER &
                                 .mint_class_rep_charge.ToString

                    srOutFileFeatures.WriteLine(strLineout)
                End With
            Next intIndex

        Catch ex As Exception
            LogErrors("FindLCMSFeatures", "Error writing to the output file (" & Path.GetFileName(strFeaturesOutputFilePath) & ")", ex, False, False, True, eLCMSFeatureFinderErrorCodes.InputFileReadError)
        Finally
            If Not srOutFileFeatures Is Nothing Then
                srOutFileFeatures.Close()
            End If
        End Try

        Try
            ' Make sure the feature to isotope map values are sorted by feature index and then by data index (this is most likely already true, but it's quick to make sure)
            SortFeatureToIsotopeMapInfo(intUMCIndex, intIsotopePeaksIndex)

            ' Write the Feature to Isotopic Peak mapping to the output file
            For intIndex = 0 To intIsotopePeaksIndex.Length - 1
                strLineout = intUMCIndex(intIndex).ToString & COLUMN_DELIMITER & intIsotopePeaksIndex(intIndex).ToString()
                srOutFileFeatureToIsotopePeakMap.WriteLine(strLineout)
            Next intIndex

        Catch ex As Exception
            LogErrors("FindLCMSFeatures", "Error writing to the output file (" & Path.GetFileName(strFeaturesOutputFilePath) & ")", ex, False, False, True, eLCMSFeatureFinderErrorCodes.InputFileReadError)
        Finally
            If Not srOutFileFeatureToIsotopePeakMap Is Nothing Then
                srOutFileFeatureToIsotopePeakMap.Close()
            End If
        End Try

        Return blnSuccess

    End Function

    Public Function FindLCMSFeaturesFromMemory(
      <Out()> ByRef objFeatures() As clsUMC,
      <Out()> ByRef intIsotopePeaksIndex() As Integer,
      <Out()> ByRef intUMCIndex() As Integer) As Boolean

        ' Populate mIsotopePeaks prior to calling this function, either by calling FindLCMSFeaturesFromFile
        '  or using AddIsotopePeak()
        ' This function will set the UMCCreator options using mFeatureFindingOptions

        Dim intMappingCount As Integer
        Dim blnSuccess As Boolean

        Try
            ' Assume success for now
            blnSuccess = True

            If mIsotopePeaks Is Nothing OrElse mIsotopePeaksCount = 0 Then
                If mIsotopePeaks Is Nothing Then
                    LogErrors("FindLCMSFeaturesFromMemory", "mIsotopePeaks is not initialized; unable to continue", Nothing, True, False, False)
                Else
                    LogErrors("FindLCMSFeaturesFromMemory", "mIsotopePeaks is empty; unable to continue", Nothing, True, False, False)
                End If

                ReDim objFeatures(-1)
                ReDim intIsotopePeaksIndex(-1)
                ReDim intUMCIndex(-1)
                Return False
            Else
                ' Shrink mIsotopePeaks to length mIsotopePeaksCount if needed
                If mIsotopePeaksCount < mIsotopePeaks.Length Then
                    ReDim Preserve mIsotopePeaks(mIsotopePeaksCount - 1)
                    ReDim Preserve mIsotopePeaksMapToIndexInDataSource(mIsotopePeaksCount - 1)
                End If
            End If

            ' Use clsUMCCreator to find LC-MS features from the loaded data
            If mUMCCreator Is Nothing Then
                mUMCCreator = New clsUMCCreator()
            Else
                mUMCCreator.ResetStatus()
            End If

            ' Set the options
            With mFeatureFindingOptions
                mUMCCreator.SetOptionsEx(
                                .MonoMassWeight,
                                .MonoMassConstraint,
                                .MonoMassConstraintIsPPM,
                                .AvgMassWeight,
                                .AvgMassConstraint,
                                .AvgMassConstraintIsPPM,
                                .LogAbundanceWeight,
                                .ScanWeight,
                                .NETWeight,
                                .FitWeight,
                                .MaxDistance,
                                .UseGenericNET,
                                .IMSDriftTimeWeight,
                                .RequireMatchingChargeState)

                mUMCCreator.MinUMCLength = .MinFeatureLengthPoints

                mUMCCreator.SetIsotopePeaks(mIsotopePeaks)

                If (.MinScan > 0 Or .MaxScan > 0) AndAlso .MinScan <> .MaxScan Then
                    mUMCCreator.SetLCMinMaxScans(.MinScan, .MaxScan)

                Else
                    Dim intMinScan = mUMCCreator.MinScan
                    Dim intMaxScan = mUMCCreator.MaxScan
                End If
            End With

            Dim objThread = New Thread(AddressOf FindLCMSFeaturesWork)

            objThread.Start()

            Do While objThread.ThreadState = ThreadState.Running Or objThread.ThreadState = ThreadState.Unstarted
                Thread.Sleep(100)
                Console.WriteLine(mUMCCreator.PercentComplete & ": " & mUMCCreator.Message)
                If mAbortProcessing Then
                    Exit Do
                End If
            Loop

            If mAbortProcessing Then
                LogErrors("FindLCMSFeaturesFromMemory", "Processing Aborted", Nothing, False, False, True, eLCMSFeatureFinderErrorCodes.UnspecifiedError)
                Try
                    objThread.Abort()
                Catch ex As Exception
                    ' Ignore errors here
                End Try
                blnSuccess = False
            Else
                blnSuccess = True
            End If

        Catch ex As Exception
            LogErrors("FindLCMSFeaturesFromMemory", "Error finding the LC-MS features using UMCCreation.dll", ex, False, False, True, eLCMSFeatureFinderErrorCodes.LCMSProcessingError)
            blnSuccess = False
        End Try

        If Not blnSuccess Then
            ReDim objFeatures(-1)
            ReDim intIsotopePeaksIndex(-1)
            ReDim intUMCIndex(-1)
            Return False
        End If

        Try
            intMappingCount = mUMCCreator.GetUmcMapping(intIsotopePeaksIndex, intUMCIndex)
            objFeatures = mUMCCreator.GetUMCs()

            ' Update the values in intIsotopePeaksIndex using mIsotopePeaksMapToIndexInDataSource
            For intIndex As Integer = 0 To intMappingCount - 1
                If intIsotopePeaksIndex(intIndex) >= 0 And intIsotopePeaksIndex(intIndex) < mIsotopePeaksMapToIndexInDataSource.Length Then
                    If intIsotopePeaksIndex(intIndex) <> mIsotopePeaksMapToIndexInDataSource(intIsotopePeaksIndex(intIndex)) Then
                        ' This code will only be reached if the source data did not start at index 0 increment up via 1 for each data point
                        intIsotopePeaksIndex(intIndex) = mIsotopePeaksMapToIndexInDataSource(intIsotopePeaksIndex(intIndex))
                    End If
                End If
            Next

            If intIsotopePeaksIndex Is Nothing Then
                ReDim intIsotopePeaksIndex(-1)
                ReDim intUMCIndex(-1)
            Else
                If intIsotopePeaksIndex.Length > intMappingCount OrElse
                   intUMCIndex.Length > intMappingCount Then
                    ReDim Preserve intIsotopePeaksIndex(intMappingCount - 1)
                    ReDim Preserve intUMCIndex(intMappingCount - 1)
                End If
            End If

            If objFeatures Is Nothing Then
                ReDim objFeatures(-1)
            End If

        Catch ex As Exception
            LogErrors("FindLCMSFeaturesFromMemory", "Error finding the LC-MS features using UMCCreation.dll", ex, False, False, True, eLCMSFeatureFinderErrorCodes.LCMSProcessingError)
            ReDim objFeatures(-1)
            ReDim intIsotopePeaksIndex(-1)
            ReDim intUMCIndex(-1)
            Return False
        End Try

        Return blnSuccess

    End Function

    Private Sub FindLCMSFeaturesWork()
        mUMCCreator.FindUMCs()
    End Sub

    Private Function GenerateOutputFileNameForLCMSFeatures(strInputFilePath As String) As String

        Dim strInputFileName As String = "UnknownDataset"

        Try
            If Not strInputFilePath Is Nothing AndAlso strInputFilePath.Length > 0 Then
                strInputFileName = Path.GetFileNameWithoutExtension(strInputFilePath)
            End If
        Catch ex As Exception
            ' Ignore errors here
        End Try

        Return strInputFileName & FILE_SUFFIX_LCMS_FEATURES & ".txt"

    End Function

    Private Function GenerateOutputFileNameForPeakToFeatureMap(strLCMSFeaturesFilePath As String) As String

        Dim strInputFileName As String = "UnknownDataset"

        Try
            If Not strLCMSFeaturesFilePath Is Nothing AndAlso strLCMSFeaturesFilePath.Length > 0 Then
                strInputFileName = Path.GetFileNameWithoutExtension(strLCMSFeaturesFilePath)

                If strInputFileName.ToLower.EndsWith(FILE_SUFFIX_LCMS_FEATURES.ToLower) Then
                    strInputFileName = strInputFileName.Substring(0, strInputFileName.Length - FILE_SUFFIX_LCMS_FEATURES.Length)
                End If
            End If
        Catch ex As Exception
            ' Ignore errors here
        End Try

        Return strInputFileName & FILE_SUFFIX_PEAK_TO_FEATURE_MAP & ".txt"


    End Function

    Private Function GetAppFolderPath() As String
        ' Could use Application.StartupPath, but .GetExecutingAssembly is better
        Return Path.GetDirectoryName(Assembly.GetExecutingAssembly().Location)
    End Function

    Private Function GetColumnValueDbl(ByRef strData() As String, intColumnIndex As Integer, Optional ByVal dblDefaultValue As Double = 0) As Double
        Try
            If intColumnIndex >= 0 Then
                Return CDbl(strData(intColumnIndex))
            Else
                Return dblDefaultValue
            End If
        Catch ex As Exception
            Return 0
        End Try
    End Function

    Private Function GetColumnValueInt(ByRef strData() As String, intColumnIndex As Integer, Optional ByVal intDefaultValue As Integer = 0) As Integer
        Try
            If intColumnIndex >= 0 Then
                Return CInt(strData(intColumnIndex))
            Else
                Return intDefaultValue
            End If
        Catch ex As Exception
            Return 0
        End Try
    End Function

    Private Function GetColumnValueSng(ByRef strData() As String, intColumnIndex As Integer, Optional ByVal sngDefaultValue As Single = 0) As Single
        Try
            If intColumnIndex >= 0 Then
                Return CSng(strData(intColumnIndex))
            Else
                Return sngDefaultValue
            End If
        Catch ex As Exception
            Return 0
        End Try
    End Function

    Public Function GetErrorMessage() As String
        ' Returns String.Empty if no error

        Dim strErrorMessage As String

        Select Case mErrorCode
            Case eLCMSFeatureFinderErrorCodes.NoError
                strErrorMessage = String.Empty
            Case eLCMSFeatureFinderErrorCodes.InvalidInputFilePath
                strErrorMessage = "Invalid input file path"
            Case eLCMSFeatureFinderErrorCodes.InvalidOutputFolderPath
                strErrorMessage = "Invalid output folder path"
            Case eLCMSFeatureFinderErrorCodes.ParameterFileNotFound
                strErrorMessage = "Parameter file not found"
            Case eLCMSFeatureFinderErrorCodes.FilePathError
                strErrorMessage = "General file path error"

            Case eLCMSFeatureFinderErrorCodes.ParameterFileReadError
                strErrorMessage = "Parameter file read error"
            Case eLCMSFeatureFinderErrorCodes.UnknownFileExtension
                strErrorMessage = "Unknown file extension"
            Case eLCMSFeatureFinderErrorCodes.InputFileReadError
                strErrorMessage = "Input file read error"
            Case eLCMSFeatureFinderErrorCodes.OutputFileWriteError
                strErrorMessage = "Error writing output file"

            Case eLCMSFeatureFinderErrorCodes.LCMSProcessingError
                strErrorMessage = "Error finding the LC-MS features using UMCCreation.dll"

            Case eLCMSFeatureFinderErrorCodes.UnspecifiedError
                strErrorMessage = "Unspecified localized error"
            Case Else
                ' This shouldn't happen
                strErrorMessage = "Unknown error state"
        End Select

        Return strErrorMessage

    End Function

    Protected Function GetParentFolderPathForFile(strFilePath As String) As String
        Dim ioFile As FileInfo

        ioFile = New FileInfo(strFilePath)
        If ioFile.Directory.Exists Then
            Return ioFile.DirectoryName
        Else
            ' Use the current working directory
            Return Path.GetDirectoryName(Assembly.GetExecutingAssembly().Location)
        End If

    End Function

    Private Function GetTempFolderPath(blnUseSystemTempPath As Boolean, blnCreateFolder As Boolean) As String
        ' If blnUseSystemTempPath = True then uses the system-defined temporary folder path
        ' Otherwise, uses the application folder as the base folder
        ' If blnCreateFolder = True, then creates the folder; if creation fails, then tries the other available folder (depending on blnUseSystemTempPath)

        Dim strFolderName As String
        Dim strFolderPathBase As String
        Dim strFolderPath As String = String.Empty
        Dim strFolderPathTest As String

        Dim intLoop As Integer
        Dim blnSuccess As Boolean

        strFolderName = Guid.NewGuid.ToString
        strFolderName = strFolderName.Replace("-"c, "_"c)

        For intLoop = 1 To 2
            If blnUseSystemTempPath Then
                strFolderPathBase = Path.GetTempPath()
            Else
                strFolderPathBase = GetAppFolderPath()
            End If

            strFolderPathTest = Path.Combine(strFolderPathBase, strFolderName)
            If intLoop = 1 Then
                strFolderPath = String.Copy(strFolderPathTest)
            End If

            blnSuccess = True
            If blnCreateFolder Then
                Try
                    Directory.CreateDirectory(strFolderPathTest)

                    ' Folder successfully created
                    strFolderPath = String.Copy(strFolderPathTest)

                Catch ex As Exception
                    ' Folder creation failed; set blnSuccess to False so that we loop
                    blnSuccess = False
                End Try
            End If

            If blnSuccess Then
                Exit For
            Else
                blnUseSystemTempPath = Not blnUseSystemTempPath
            End If
        Next intLoop

        Return strFolderPath
    End Function

    Private Sub DefineDefaultColumnMappings(ByRef intColumnMapping() As Integer)

        ReDim intColumnMapping(INPUT_FILE_COLUMN_NAME_COUNT - 1)

        intColumnMapping(eInputFileColumnNames.Scan) = eInputFileColumnNames.Scan
        intColumnMapping(eInputFileColumnNames.Charge) = eInputFileColumnNames.Charge
        intColumnMapping(eInputFileColumnNames.Abundance) = eInputFileColumnNames.Abundance
        intColumnMapping(eInputFileColumnNames.MZOfMostAbu) = eInputFileColumnNames.MZOfMostAbu
        intColumnMapping(eInputFileColumnNames.Fit) = eInputFileColumnNames.Fit
        intColumnMapping(eInputFileColumnNames.AverageMass) = eInputFileColumnNames.AverageMass
        intColumnMapping(eInputFileColumnNames.MonoisotopicMass) = eInputFileColumnNames.MonoisotopicMass
        intColumnMapping(eInputFileColumnNames.MassOfMostAbundant) = eInputFileColumnNames.MassOfMostAbundant
        intColumnMapping(eInputFileColumnNames.FWHM) = eInputFileColumnNames.FWHM
        intColumnMapping(eInputFileColumnNames.SignalToNoise) = eInputFileColumnNames.SignalToNoise
        intColumnMapping(eInputFileColumnNames.MonoMassAbundance) = eInputFileColumnNames.MonoMassAbundance
        intColumnMapping(eInputFileColumnNames.MonoMassPlus2DaAbundance) = eInputFileColumnNames.MonoMassPlus2DaAbundance
        intColumnMapping(eInputFileColumnNames.Index) = eInputFileColumnNames.Index
        intColumnMapping(eInputFileColumnNames.NET) = eInputFileColumnNames.NET
        intColumnMapping(eInputFileColumnNames.IMSDriftTime) = eInputFileColumnNames.IMSDriftTime

    End Sub

    Private Sub InitializeColumnMappings(<Out()> ByRef intColumnMapping() As Integer)
        Dim intIndex As Integer

        ReDim intColumnMapping(INPUT_FILE_COLUMN_NAME_COUNT - 1)
        For intIndex = 0 To intColumnMapping.Length - 1
            intColumnMapping(intIndex) = -1
        Next

    End Sub

    Private Sub InitializeLocalVariables()
        mErrorCode = eLCMSFeatureFinderErrorCodes.NoError
        mStatusMessage = String.Empty

        mProgressStepDescription = String.Empty
        mProgressPercentComplete = 0

        SetDefaultFeatureFindingOptions(mFeatureFindingOptions)

        ClearIsotopePeaks()
    End Sub

    Public Shared Function IsNumber(strValue As String) As Boolean
        Dim result As Double
        Try
            Return Double.TryParse(strValue, result)
        Catch ex As Exception
            Return False
        End Try
    End Function

    ' This function is not currently implemented
    ''
    ''Public Function LoadParameterFileSettings(ByVal strParameterFilePath As String) As Boolean

    ''    Dim objSettingsFile As New PRISM.Files.XmlSettingsFileAccessor

    ''    Try

    ''        If strParameterFilePath Is Nothing OrElse strParameterFilePath.Length = 0 Then
    ''            ' No parameter file specified; nothing to load
    ''            Return True
    ''        End If

    ''        If Not System.IO.File.Exists(strParameterFilePath) Then
    ''            ' See if strParameterFilePath points to a file in the same directory as the application
    ''            strParameterFilePath = System.IO.Path.Combine(GetAppFolderPath(), System.IO.Path.GetFileName(strParameterFilePath))
    ''            If Not System.IO.File.Exists(strParameterFilePath) Then
    ''                LogErrors("LoadParameterFileSettings", "Parameter file not found: " & strParameterFilePath, Nothing, True, False, False)
    ''                SetErrorCode(eLCMSFeatureFinderErrorCodes.ParameterFileNotFound)
    ''                Return False
    ''            End If
    ''        End If

    ''        ' Pass False to .LoadSettings() here to turn off case sensitive matching
    ''        If objSettingsFile.LoadSettings(strParameterFilePath, False) Then
    ''            With objSettingsFile

    ''                If Not .SectionPresent(XML_SECTION_LCMS_FEATURE_FINDER_SETTINGS) Then
    ''                    ' MS File Scanner section not found; that's ok
    ''                Else
    ''                    Me.IncludeAssociatedMS1Spectra = .GetParam(XML_SECTION_LCMS_FEATURE_FINDER_SETTINGS, "FutureAppSetting", Me.FutureAppSetting)
    ''                End If

    ''            End With
    ''        Else
    ''            LogErrors("LoadParameterFileSettings", "Error calling objSettingsFile.LoadSettings for " & strParameterFilePath, Nothing, True, False, True, eLCMSFeatureFinderErrorCodes.ParameterFileReadError)
    ''            Return False
    ''        End If

    ''    Catch ex As System.Exception
    ''        LogErrors("LoadParameterFileSettings", "Error in LoadParameterFileSettings", ex, True, False, True, eLCMSFeatureFinderErrorCodes.ParameterFileReadError)
    ''        Return False
    ''    End Try

    ''    Return True

    ''End Function

    Private Sub LogErrors(strSource As String, strMessage As String, ex As Exception, Optional ByVal blnAllowInformUser As Boolean = True, Optional ByVal blnAllowThrowingException As Boolean = True, Optional ByVal blnLogLocalOnly As Boolean = True, Optional ByVal eNewErrorCode As eLCMSFeatureFinderErrorCodes = eLCMSFeatureFinderErrorCodes.NoError)
        Dim strMessageWithoutCRLF As String

        mStatusMessage = String.Copy(strMessage)

        strMessageWithoutCRLF = mStatusMessage.Replace(ControlChars.NewLine, "; ")

        If ex Is Nothing Then
            ex = New Exception("Error")
        Else
            If Not ex.Message Is Nothing AndAlso ex.Message.Length > 0 Then
                strMessageWithoutCRLF &= "; " & ex.Message
            End If
        End If

        Console.WriteLine(DateTime.Now().ToLongTimeString & "; " & strMessageWithoutCRLF, strSource)

        'If Not mErrorLogger Is Nothing Then
        '    mErrorLogger.PostError(mStatusMessage.Replace(ControlChars.NewLine, "; "), ex, blnLogLocalOnly)
        'End If

        If Not eNewErrorCode = eLCMSFeatureFinderErrorCodes.NoError Then
            SetErrorCode(eNewErrorCode, True)
        End If

        If blnAllowThrowingException Then
            Throw New Exception(mStatusMessage, ex)
        End If
    End Sub

    Protected Sub OperationComplete()
        RaiseEvent ProgressComplete()
    End Sub

    Private Sub ParseInputFileColumnNames(strSplitLine() As String, ByRef intColumnMapping() As Integer)
        Dim intIndex As Integer
        Dim strCurrentColumn As String

        Dim htColumnInfo As Hashtable
        Dim objEnum As IDictionaryEnumerator

        Try
            htColumnInfo = New Hashtable

            For intIndex = 0 To strSplitLine.Length - 1
                strCurrentColumn = strSplitLine(intIndex).Trim
                If Not htColumnInfo.Contains(strCurrentColumn) Then
                    htColumnInfo.Add(strCurrentColumn, intIndex)
                End If
            Next intIndex

            If htColumnInfo.Count > 0 Then
                If intColumnMapping Is Nothing Then
                    ReDim intColumnMapping(INPUT_FILE_COLUMN_NAME_COUNT - 1)
                End If

                objEnum = htColumnInfo.GetEnumerator
                Do While objEnum.MoveNext
                    strCurrentColumn = CStr(objEnum.Key).ToLower

                    Select Case strCurrentColumn
                        Case ISOS_COLUMN_SCAN_NUM
                            intColumnMapping(eInputFileColumnNames.Scan) = CInt(objEnum.Value)
                        Case ISOS_COLUMN_FRAME_NUM
                            intColumnMapping(eInputFileColumnNames.Scan) = CInt(objEnum.Value)
                        Case ISOS_COLUMN_IMS_SCAN_NUM
                            ' This is a valid column, but LC-MS Feature finder doesn't use it
                        Case ISOS_COLUMN_CHARGE
                            intColumnMapping(eInputFileColumnNames.Charge) = CInt(objEnum.Value)
                        Case ISOS_COLUMN_ABUNDANCE
                            intColumnMapping(eInputFileColumnNames.Abundance) = CInt(objEnum.Value)
                        Case ISOS_COLUMN_MZ
                            intColumnMapping(eInputFileColumnNames.MZOfMostAbu) = CInt(objEnum.Value)
                        Case ISOS_COLUMN_FIT
                            intColumnMapping(eInputFileColumnNames.Fit) = CInt(objEnum.Value)
                        Case ISOS_COLUMN_AVERAGE_MW
                            intColumnMapping(eInputFileColumnNames.AverageMass) = CInt(objEnum.Value)
                        Case ISOS_COLUMN_MONOISOTOPIC_MW
                            intColumnMapping(eInputFileColumnNames.MonoisotopicMass) = CInt(objEnum.Value)
                        Case ISOS_COLUMN_MOSTABUNDANT_MW
                            intColumnMapping(eInputFileColumnNames.MassOfMostAbundant) = CInt(objEnum.Value)
                        Case ISOS_COLUMN_FWHM
                            intColumnMapping(eInputFileColumnNames.FWHM) = CInt(objEnum.Value)
                        Case ISOS_COLUMN_SIGNAL_NOISE
                            intColumnMapping(eInputFileColumnNames.SignalToNoise) = CInt(objEnum.Value)
                        Case ISOS_COLUMN_MONO_ABUNDANCE
                            intColumnMapping(eInputFileColumnNames.MonoMassAbundance) = CInt(objEnum.Value)
                        Case ISOS_COLUMN_MONO_PLUS2_ABUNDANCE
                            intColumnMapping(eInputFileColumnNames.MonoMassPlus2DaAbundance) = CInt(objEnum.Value)
                        Case ISOS_COLUMN_INDEX
                            intColumnMapping(eInputFileColumnNames.Index) = CInt(objEnum.Value)
                        Case ISOS_COLUMN_NET
                            intColumnMapping(eInputFileColumnNames.NET) = CInt(objEnum.Value)
                        Case ISOS_COLUMN_IMS_DRIFT_TIME
                            intColumnMapping(eInputFileColumnNames.IMSDriftTime) = CInt(objEnum.Value)
                        Case ISOS_COLUMN_IMS_DRIFT_TIME_ALT
                            intColumnMapping(eInputFileColumnNames.IMSDriftTime) = CInt(objEnum.Value)
                        Case ISOS_COLUMN_ORIG_INTENSITY
                            ' This is a valid column, but LC-MS Feature finder doesn't use it
                        Case ISOS_COLUMN_TIA_ORIG_INTENSITY
                            ' This is a valid column, but LC-MS Feature finder doesn't use it
                        Case Else
                            ' Ignore this column (display a warning on the console)
                            Console.WriteLine("Ignoring column '" & strCurrentColumn & "' since the name is not recognized")
                    End Select
                Loop
            End If
        Catch ex As Exception
            ' Error parsing column names; leave
            LogErrors("ParseInputFileColumnNames", "Warning: Error parsing the column names on the first line of the input file", ex, True, False, False, eLCMSFeatureFinderErrorCodes.InputFileReadError)
        End Try

    End Sub

    Public Function ProcessFile(strInputFilePath As String, strFeaturesOutputFilePath As String, blnResetErrorCode As Boolean) As Boolean
        ' strInputFilePath should specify a text file with LC-MS data
        ' strFeaturesOutputFilePath should specify the file in which to write the LC-MS features; note that the FeatureToPeakMapFilePath will be auto-defined based on strFeaturesOutputFilePath
        ' strFeaturesOutputFilePath will be auto-defined if blank
        '
        ' This function returns True if successful, false if an error

        ' This function will look for a .Ini file with the same name as strInputFilePath but ending in .ini
        ' If found, then opens it to read the UMC Creation settings

        Dim strInputFolderPath As String
        Dim strOutputFolderPath As String

        Dim strFeatureToPeakMapFilePath As String

        Dim blnSuccess As Boolean

        Try

            If strInputFilePath Is Nothing OrElse strInputFilePath.Length = 0 Then
                ' Nothing to do
                LogErrors("ProcessFile", "Input filename not defined; unable to continue", Nothing, True, False, True, eLCMSFeatureFinderErrorCodes.InvalidInputFilePath)
                Return False
            End If

            strInputFolderPath = GetParentFolderPathForFile(strInputFilePath)

            If strFeaturesOutputFilePath Is Nothing OrElse strFeaturesOutputFilePath.Length = 0 Then
                ' Auto-define strFeaturesOutputFilePath
                strFeaturesOutputFilePath = Path.Combine(strInputFolderPath, GenerateOutputFileNameForLCMSFeatures(strInputFilePath))
            End If

            ' Always auto-define strFeatureToPeakMapFilePath
            strOutputFolderPath = GetParentFolderPathForFile(strFeaturesOutputFilePath)
            strFeatureToPeakMapFilePath = Path.Combine(strOutputFolderPath, GenerateOutputFileNameForPeakToFeatureMap(strFeaturesOutputFilePath))

        Catch ex As Exception
            LogErrors("ProcessFile", "Error validating the filenames", ex, True, False, False, eLCMSFeatureFinderErrorCodes.FilePathError)
            Return False
        End Try

        Try
            If Not File.Exists(strInputFilePath) Then
                LogErrors("ProcessFile", "Input file not found: " & strInputFilePath, Nothing, True, False, True, eLCMSFeatureFinderErrorCodes.InvalidInputFilePath)
                Return False
            End If

            ' Look for a .Ini file matching strInputFilePath
            ' If found, load the options
            AutoLoadOptions(strInputFilePath)

            blnSuccess = FindLCMSFeaturesFromFile(strInputFilePath, strFeaturesOutputFilePath, strFeatureToPeakMapFilePath)

        Catch ex As Exception
            LogErrors("ProcessFile", "Error calling ExtractScanStatsFromXcaliburDataFile", ex, True, False, False, eLCMSFeatureFinderErrorCodes.UnspecifiedError)
            Return False
        End Try

        Return blnSuccess

    End Function

    Public Function ProcessFileWildcard(strInputFilePath As String, strOutputFileOrFolderPath As String, blnResetErrorCode As Boolean) As Boolean
        ' Returns True if success, False if failure

        Dim blnSuccess As Boolean
        Dim intMatchCount As Integer

        Dim strCleanPath As String
        Dim strInputFolderPath As String
        Dim strOutputFilePath As String
        Dim strOutputFolderPath As String

        Dim ioFileMatch As FileInfo

        Dim ioFileInfo As FileInfo
        Dim ioInputFolderInfo As DirectoryInfo
        Dim ioOutputFolderInfo As DirectoryInfo = Nothing


        Dim intIndex As Integer
        Dim strFilePaths() As String

        mAbortProcessing = False
        blnSuccess = True
        Try
            ' Possibly reset the error code
            If blnResetErrorCode Then
                SetErrorCode(eLCMSFeatureFinderErrorCodes.NoError)
            End If

            ' See if strOutputFileOrFolderPath is a folder
            Try
                If Not strOutputFileOrFolderPath Is Nothing AndAlso strOutputFileOrFolderPath.Length > 0 Then
                    If Directory.Exists(strOutputFileOrFolderPath) Then
                        ioOutputFolderInfo = New DirectoryInfo(strOutputFileOrFolderPath)
                    End If
                End If
            Catch ex As Exception
                ioOutputFolderInfo = Nothing
            End Try

            ' See if strInputFilePath contains a wildcard
            If Not strInputFilePath Is Nothing AndAlso (strInputFilePath.IndexOf("*"c) >= 0 Or strInputFilePath.IndexOf("?"c) >= 0) Then
                ' Obtain a list of the matching files and folders

                ' Copy the path into strCleanPath and replace any * or ? characters with
                strCleanPath = strInputFilePath.Replace("*", "_")
                strCleanPath = strCleanPath.Replace("?", "_")

                ioFileInfo = New FileInfo(strCleanPath)
                If ioFileInfo.Directory.Exists Then
                    strInputFolderPath = ioFileInfo.DirectoryName
                Else
                    ' Use the current working directory
                    strInputFolderPath = Path.GetDirectoryName(Assembly.GetExecutingAssembly().Location)
                End If

                ioInputFolderInfo = New DirectoryInfo(strInputFolderPath)

                ' Remove any directory information from strinputFilePath
                strInputFilePath = Path.GetFileName(strInputFilePath)

                ' Cache the list of input files now before starting

                intMatchCount = 0
                ReDim strFilePaths(4)

                intMatchCount = 0
                For Each ioFileMatch In ioInputFolderInfo.GetFiles(strInputFilePath)
                    If intMatchCount >= strFilePaths.Length Then
                        ReDim Preserve strFilePaths(strFilePaths.Length * 2 - 1)
                    End If

                    strFilePaths(intMatchCount) = ioFileMatch.FullName
                    intMatchCount += 1
                Next ioFileMatch

                ' Now process each of the input files
                For intIndex = 0 To intMatchCount - 1
                    strOutputFilePath = GenerateOutputFileNameForLCMSFeatures(strFilePaths(intIndex))

                    If ioOutputFolderInfo Is Nothing Then
                        strOutputFilePath = Path.Combine(ioInputFolderInfo.FullName, strOutputFilePath)
                    Else
                        strOutputFilePath = Path.Combine(ioOutputFolderInfo.FullName, strOutputFilePath)
                    End If

                    blnSuccess = ProcessFile(strFilePaths(intIndex), strOutputFilePath, blnResetErrorCode)

                    If Not blnSuccess Or mAbortProcessing Then Exit For
                Next intIndex

                If intMatchCount = 0 Then
                    If mErrorCode = eLCMSFeatureFinderErrorCodes.NoError Then
                        LogErrors("ProcessFileWildcard", "No match was found for the input file path:" & strInputFilePath, Nothing, True, False, False, eLCMSFeatureFinderErrorCodes.InvalidInputFilePath)
                    End If
                End If
            Else
                If ioOutputFolderInfo Is Nothing Then
                    strOutputFolderPath = String.Copy(strOutputFileOrFolderPath)
                Else
                    strOutputFolderPath = Path.Combine(ioOutputFolderInfo.FullName, GenerateOutputFileNameForLCMSFeatures(strInputFilePath))
                End If

                blnSuccess = ProcessFile(strInputFilePath, strOutputFolderPath, blnResetErrorCode)
            End If

        Catch ex As Exception
            LogErrors("ProcessFileWildcard", "Error processing files:", ex, True, True, False, eLCMSFeatureFinderErrorCodes.InvalidInputFilePath)
        End Try

        Return blnSuccess

    End Function

    Protected Sub ResetProgress()
        RaiseEvent ProgressReset()
    End Sub

    Protected Sub ResetProgress(strProgressStepDescription As String)
        UpdateProgress(strProgressStepDescription, 0)
        RaiseEvent ProgressReset()
    End Sub

    ' This function is not currently implemented
    ''
    ''Public Function SaveParameterFileSettings(ByVal strParameterFilePath As String) As Boolean

    ''    Dim objSettingsFile As New PRISM.Files.XmlSettingsFileAccessor

    ''    Dim intIndex As Integer

    ''    Try

    ''        If strParameterFilePath Is Nothing OrElse strParameterFilePath.Length = 0 Then
    ''            ' No parameter file specified; unable to save
    ''            Return False
    ''        End If

    ''        ' Pass True to .LoadSettings() here so that newly made Xml files will have the correct capitalization
    ''        If objSettingsFile.LoadSettings(strParameterFilePath, True) Then
    ''            With objSettingsFile

    ''                ' General settings
    ''                .SetParam(XML_SECTION_LCMS_FEATURE_FINDER_SETTINGS, "FutureAppSetting", Me.FutureAppSetting)
    ''            End With

    ''            objSettingsFile.SaveSettings()

    ''        End If

    ''    Catch ex As System.Exception
    ''        LogErrors("SaveParameterFileSettings", "Error in SaveParameterFileSettings", ex, True, False, False, eLCMSFeatureFinderErrorCodes.OutputFileWriteError)
    ''        Return False
    ''    Finally
    ''        objSettingsFile = Nothing
    ''    End Try

    ''    Return True

    ''End Function

    Private Sub SetDefaultFeatureFindingOptions(ByRef udtFeatureFindingOptions As udtFeatureFindingOptionsType)
        With udtFeatureFindingOptions
            .MonoMassWeight = 0.01
            .MonoMassConstraint = 10
            .MonoMassConstraintIsPPM = True

            .AvgMassWeight = 0.01
            .AvgMassConstraint = 10
            .AvgMassConstraintIsPPM = True

            .LogAbundanceWeight = 0.1
            .ScanWeight = 0.005         ' Used if .UseGenericNET = False
            .NETWeight = 15             ' Used if .UseGenericNET = True
            .FitWeight = 0.1
            .IMSDriftTimeWeight = 0.1

            .MaxDistance = 0.1
            .UseGenericNET = True

            .MinScan = 0
            .MaxScan = 0

            .MinFeatureLengthPoints = 2
            .RequireMatchingChargeState = False
        End With
    End Sub

    Private Sub SetErrorCode(eNewErrorCode As eLCMSFeatureFinderErrorCodes)
        SetErrorCode(eNewErrorCode, False)
    End Sub

    Private Sub SetErrorCode(eNewErrorCode As eLCMSFeatureFinderErrorCodes, blnLeaveExistingErrorCodeUnchanged As Boolean)

        If blnLeaveExistingErrorCodeUnchanged AndAlso mErrorCode <> eLCMSFeatureFinderErrorCodes.NoError Then
            ' An error code is already defined; do not change it
        Else
            mErrorCode = eNewErrorCode
        End If

    End Sub

    Private Sub SortFeatureToIsotopeMapInfo(ByRef intUMCIndex() As Integer, ByRef intIsotopePeaksIndex() As Integer)

        Dim intCurrentUMC As Integer

        Dim intStartIndex As Integer
        Dim intLastStartIndex As Integer

        Dim intEndIndex As Integer
        Dim intLastIndex As Integer

        ' Sort intUMCIndex() and intIsotopePeaksIndex() in parallel
        Array.Sort(intUMCIndex, intIsotopePeaksIndex)

        ' Fine-tune the sort such that the peaks for each UMC are sorted ascending by index
        intStartIndex = 0
        intLastStartIndex = -1
        intEndIndex = 0

        Try
            intLastIndex = intUMCIndex.Length - 1
            Do While intStartIndex < intLastIndex
                intCurrentUMC = intUMCIndex(intStartIndex)

                If intLastStartIndex = intStartIndex Then
                    ' This shouldn't normally happen
                    Exit Do
                Else
                    intLastStartIndex = intStartIndex
                End If

                Do While intEndIndex < intLastIndex AndAlso intUMCIndex(intEndIndex + 1) = intCurrentUMC
                    intEndIndex += 1
                Loop

                If intEndIndex > intStartIndex Then
                    ' Sort the data for this UMC
                    Array.Sort(intIsotopePeaksIndex, intStartIndex, intEndIndex - intStartIndex + 1)
                End If

                intStartIndex = intEndIndex + 1
            Loop

        Catch ex As Exception
            Console.WriteLine("Error in SortFeatureToIsotopeMapInfo: " & ex.Message)
        End Try

    End Sub

    Protected Sub UpdateProgress(strProgressStepDescription As String)
        UpdateProgress(strProgressStepDescription, mProgressPercentComplete)
    End Sub

    Protected Sub UpdateProgress(dblPercentComplete As Double)
        UpdateProgress(Me.ProgressStepDescription, CSng(dblPercentComplete))
    End Sub

    Protected Sub UpdateProgress(sngPercentComplete As Single)
        UpdateProgress(Me.ProgressStepDescription, sngPercentComplete)
    End Sub

    Protected Sub UpdateProgress(strProgressStepDescription As String, sngPercentComplete As Single)
        mProgressStepDescription = String.Copy(strProgressStepDescription)
        If sngPercentComplete < 0 Then
            sngPercentComplete = 0
        ElseIf sngPercentComplete > 100 Then
            sngPercentComplete = 100
        End If
        mProgressPercentComplete = sngPercentComplete

        RaiseEvent ProgressChanged(Me.ProgressStepDescription, Me.ProgressPercentComplete)
    End Sub

    Private Function ValueToString(dblValue As Double, intDigitsOfPrecision As Integer, Optional ByVal sngScientificNotationThreshold As Single = 1000000) As String
        Dim strFormatString As String
        Dim strValue As String
        Dim strMantissa As String

        Dim dblNewValue As Double

        Dim intDigitsAfterDecimal As Integer

        If intDigitsOfPrecision < 1 Then intDigitsOfPrecision = 1

        Try
            strMantissa = "0." & New String("0"c, Math.Max(intDigitsOfPrecision - 1, 1)) & "E+00"

            If dblValue = 0 Then
                strValue = "0"
            ElseIf Math.Abs(dblValue) < 1 Then
                strFormatString = "0." & New String("0"c, intDigitsOfPrecision)
                strValue = dblValue.ToString(strFormatString)
                dblNewValue = Double.Parse(strValue)

                If dblNewValue = 0 Then
                    strValue = dblValue.ToString(strMantissa)
                Else
                    strValue = strValue.TrimEnd("0"c)
                End If
            Else
                intDigitsAfterDecimal = intDigitsOfPrecision - CInt(Math.Ceiling(Math.Log10(Math.Abs(dblValue))))

                If dblValue >= sngScientificNotationThreshold Then
                    strValue = dblValue.ToString(strMantissa)
                Else
                    If intDigitsAfterDecimal > 0 Then
                        strValue = dblValue.ToString("0." & New String("0"c, intDigitsAfterDecimal))
                        strValue = strValue.TrimEnd("0"c)
                        strValue = strValue.TrimEnd("."c)
                    Else
                        strValue = dblValue.ToString("0")
                    End If
                End If
            End If

        Catch ex As Exception
            LogErrors("ValueToString", "Error converting " & dblValue.ToString & " to a string", ex, False, False, False, eLCMSFeatureFinderErrorCodes.UnspecifiedError)
            Return dblValue.ToString
        End Try

        Return strValue

    End Function

End Class

