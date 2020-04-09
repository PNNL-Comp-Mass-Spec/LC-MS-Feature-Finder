Option Strict On

Imports System.Collections.Generic
Imports System.IO
Imports System.Reflection
Imports System.Threading

' This module reads a text file with mass and intensity data for MS spectra
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

Module modMain

    Public Const PROGRAM_DATE As String = "July 26, 2016"

    Private mInputDataFilePath As String
    Private mOutputFileOrFolderPath As String

    Private mParameterFilePath As String            ' Optional

    Private mShowExampleIni As Boolean

    Public Function Main() As Integer
        ' Returns 0 if no error, error code if an error

        Dim objLCMSFeatureFinder As clsLCMSFeatureFinder
        Dim objParseCommandLine As New clsParseCommandLine

        Dim blnProceed As Boolean
        Dim blnSuccess As Boolean

        ' Initialize default values
        Dim intReturnCode = 0

        mInputDataFilePath = String.Empty
        mOutputFileOrFolderPath = String.Empty

        mParameterFilePath = String.Empty

        Try
            blnProceed = False
            If objParseCommandLine.ParseCommandLine Then
                If SetOptionsUsingCommandLineParameters(objParseCommandLine) Then blnProceed = True
            End If

            If mShowExampleIni Then
                ShowExampleIniFile()
                Return 0
            End If

            If Not blnProceed OrElse
               objParseCommandLine.NeedToShowHelp OrElse
               mInputDataFilePath.Length = 0 Then
                ShowProgramHelp()
                intReturnCode = -1
            Else
                objLCMSFeatureFinder = New clsLCMSFeatureFinder

                With objLCMSFeatureFinder

                    ''If Not mParameterFilePath Is Nothing AndAlso mParameterFilePath.Length > 0 Then
                    ''    .LoadParameterFileSettings(mParameterFilePath)
                    ''End If
                End With

                blnSuccess = objLCMSFeatureFinder.ProcessFileWildcard(mInputDataFilePath, mOutputFileOrFolderPath, True)

                If blnSuccess Then
                    intReturnCode = 0
                Else
                    intReturnCode = objLCMSFeatureFinder.ErrorCode
                    If intReturnCode = 0 Then
                        ShowErrorMessage("Unknown error while processing; error code: " & objLCMSFeatureFinder.ErrorCode)
                    Else
                        ShowErrorMessage("Error while processing: " & objLCMSFeatureFinder.GetErrorMessage())
                    End If

                    System.Threading.Thread.Sleep(750)
                End If

            End If

        Catch ex As Exception
            ShowErrorMessage("Error occurred in modMain->Main: " & Environment.NewLine & ex.Message)
            intReturnCode = -1
        End Try

        Return intReturnCode

    End Function

    Private Function SetOptionsUsingCommandLineParameters(objParseCommandLine As clsParseCommandLine) As Boolean
        ' Returns True if no problems; otherwise, returns false

        Dim strValue As String = String.Empty
        Dim strValidParameters = New String() {"I", "O", "P", "ShowIni"}

        Try
            ' Make sure no invalid parameters are present
            If objParseCommandLine.InvalidParametersPresent(strValidParameters) Then
                Return False
            Else

                ' Query objParseCommandLine to see if various parameters are present
                With objParseCommandLine
                    If .NonSwitchParameterCount > 0 Then
                        mInputDataFilePath = .RetrieveNonSwitchParameter(0)
                    End If

                    If .RetrieveValueForParameter("I", strValue) Then mInputDataFilePath = strValue
                    If .RetrieveValueForParameter("O", strValue) Then mOutputFileOrFolderPath = strValue

                    If .RetrieveValueForParameter("P", strValue) Then mParameterFilePath = strValue

                    If .IsParameterPresent("ShowIni") Then
                        mShowExampleIni = True
                    End If

                End With

                Return True
            End If

        Catch ex As Exception
            ShowErrorMessage("Error parsing the command line parameters: " & Environment.NewLine & ex.Message)
        End Try
        Return False

    End Function

    Private Sub ShowErrorMessage(strMessage As String)
        Dim strSeparator = "------------------------------------------------------------------------------"

        Console.WriteLine()
        Console.WriteLine(strSeparator)
        Console.WriteLine(strMessage)
        Console.WriteLine(strSeparator)
        Console.WriteLine()

        WriteToErrorStream(strMessage)
    End Sub

    Private Sub ShowErrorMessage(strTitle As String, items As List(Of String))
        Dim strSeparator = "------------------------------------------------------------------------------"
        Dim strMessage As String

        Console.WriteLine()
        Console.WriteLine(strSeparator)
        Console.WriteLine(strTitle)
        strMessage = strTitle & ":"

        For Each item As String In items
            Console.WriteLine("   " + item)
            strMessage &= " " & item
        Next
        Console.WriteLine(strSeparator)
        Console.WriteLine()

        WriteToErrorStream(strMessage)
    End Sub

    Private Sub ShowExampleIniFile()
        Console.WriteLine()
        Console.WriteLine("[UMCCreationOptions]")
        Console.WriteLine("MonoMassWeight=0.01")
        Console.WriteLine("MonoMassConstraint=10")
        Console.WriteLine("MonoMassConstraintIsPPM=True")
        Console.WriteLine("AvgMassWeight=0")
        Console.WriteLine("AvgMassConstraint=10")
        Console.WriteLine("AvgMassConstraintIsPPM=True")
        Console.WriteLine("LogAbundanceWeight=0.1")
        Console.WriteLine("NETWeight=15")
        Console.WriteLine("FitWeight=0.1")
        Console.WriteLine("ScanWeight=0.005")
        Console.WriteLine("IMSDriftTimeWeight=0")
        Console.WriteLine("MaxDistance=0.1")
        Console.WriteLine("UseGenericNET=True")
        Console.WriteLine("MinScan=0")
        Console.WriteLine("MaxScan=0")
        Console.WriteLine("MinFeatureLengthPoints=2")
        Console.WriteLine("UseCharge=False")
        Console.WriteLine()
    End Sub

    Private Sub ShowProgramHelp()

        Try
            Dim exeName = Path.GetFileName(Assembly.GetExecutingAssembly().Location)

            Console.WriteLine("This program will read a text file with mass and intensity data for MS spectra and determine the LC-MS features present using UMCCreation.dll.")
            Console.WriteLine()
            Console.WriteLine("Program syntax:" & ControlChars.NewLine & exeName)
            Console.WriteLine(" /I:InputFileName.txt [/O:OutputFileName.txt] [/ShowIni]")
            Console.WriteLine()
            Console.WriteLine("The input file name is required. If the filename contains spaces, surround it with double quotes. ")
            Console.WriteLine()
            Console.WriteLine("If the output file name is not supplied, it will be auto-generated as InputFileName_Features.txt. ")
            Console.WriteLine()
            Console.WriteLine("This program will look for a .Ini file with the same name as the input file, but with extension .Ini")
            Console.WriteLine("If found, the .Ini file is processed to load the weight values to use when clustering the data.")
            Console.WriteLine("To se an example .ini File, use " & exeName & " /ShowIni")
            Console.WriteLine()

            Console.WriteLine("Program written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA) in 2007")
            Console.WriteLine("Copyright 2007, Battelle Memorial Institute.  All Rights Reserved.")
            Console.WriteLine()
            Console.WriteLine("E-mail: matthew.monroe@pnl.gov or matt@alchemistmatt.com")
            Console.WriteLine("Website: http://omics.pnl.gov/ or http://panomics.pnnl.gov/")

            Thread.Sleep(750)

        Catch ex As Exception
            ShowErrorMessage("Error displaying the program syntax: " & ex.Message)
        End Try

    End Sub

    Private Sub WriteToErrorStream(strErrorMessage As String)
        Try
            Using swErrorStream = New StreamWriter(Console.OpenStandardError())
                swErrorStream.WriteLine(strErrorMessage)
            End Using
        Catch ex As Exception
            ' Ignore errors here
        End Try
    End Sub

End Module
