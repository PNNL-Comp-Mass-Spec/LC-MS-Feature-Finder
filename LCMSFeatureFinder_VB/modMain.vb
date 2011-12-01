Option Strict On

' This module reads a text file with mass and intensity data for MS spectra
' and determines the LC-MS features present using UMCCreation.dll
'
' -------------------------------------------------------------------------------
' Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
' Program started January 30, 2007
' Copyright 2007, Battelle Memorial Institute.  All Rights Reserved.

' E-mail: matthew.monroe@pnl.gov or matt@alchemistmatt.com
' Website: http://ncrr.pnl.gov/ or http://www.sysbio.org/resources/staff/
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

	Public Const PROGRAM_DATE As String = "November 30, 2011"

    Private mInputDataFilePath As String
    Private mOutputFileOrFolderPath As String
    Private mDatasetNum As Integer

    Private mParameterFilePath As String            ' Optional

    Private mQuietMode As Boolean

    Public Function Main() As Integer
        ' Returns 0 if no error, error code if an error

        Dim intReturnCode As Integer
        Dim objLCMSFeatureFinder As clsLCMSFeatureFinder
        Dim objParseCommandLine As New clsParseCommandLine

        Dim blnProceed As Boolean
        Dim blnSuccess As Boolean

        ' Initialize default values
        intReturnCode = 0

        mInputDataFilePath = String.Empty
        mOutputFileOrFolderPath = String.Empty
        mDatasetNum = 0

        mParameterFilePath = String.Empty

        Try
            blnProceed = False
            If objParseCommandLine.ParseCommandLine Then
                If SetOptionsUsingCommandLineParameters(objParseCommandLine) Then blnProceed = True
            End If

            If objParseCommandLine.ParameterCount = 0 OrElse Not blnProceed OrElse _
               objParseCommandLine.NeedToShowHelp OrElse _
               mInputDataFilePath.Length = 0 Then
                ShowProgramHelp()
                intReturnCode = -1
            Else
                objLCMSFeatureFinder = New clsLCMSFeatureFinder

                With objLCMSFeatureFinder
                    .ShowMessages = Not mQuietMode

                    ''If Not mParameterFilePath Is Nothing AndAlso mParameterFilePath.Length > 0 Then
                    ''    .LoadParameterFileSettings(mParameterFilePath)
                    ''End If
                End With

                blnSuccess = objLCMSFeatureFinder.ProcessFileWildcard(mInputDataFilePath, mOutputFileOrFolderPath, True)

                If blnSuccess Then
                    intReturnCode = 0
                Else
                    intReturnCode = objLCMSFeatureFinder.ErrorCode
                    If intReturnCode <> 0 AndAlso Not mQuietMode Then
                        MsgBox("Error while processing: " & objLCMSFeatureFinder.GetErrorMessage(), MsgBoxStyle.Exclamation Or MsgBoxStyle.OKOnly, "Error")
                    End If
                End If

            End If

        Catch ex As System.Exception
            If mQuietMode Then
                Throw ex
            Else
                MsgBox("Error occurred: " & ControlChars.NewLine & ex.Message, MsgBoxStyle.Exclamation Or MsgBoxStyle.OKOnly, "Error")
            End If
            intReturnCode = -1
        End Try

        Return intReturnCode

    End Function

    Private Function SetOptionsUsingCommandLineParameters(ByVal objParseCommandLine As clsParseCommandLine) As Boolean
        ' Returns True if no problems; otherwise, returns false

		Dim strValue As String = String.Empty
        Dim strValidParameters() As String = New String() {"I", "O", "D", "P", "Q"}

        Try
            ' Make sure no invalid parameters are present
            If objParseCommandLine.InvalidParametersPresent(strValidParameters) Then
                Return False
            Else

                ' Query objParseCommandLine to see if various parameters are present
                With objParseCommandLine
                    If .RetrieveValueForParameter("I", strValue) Then mInputDataFilePath = strValue
                    If .RetrieveValueForParameter("O", strValue) Then mOutputFileOrFolderPath = strValue
                    If .RetrieveValueForParameter("D", strValue) Then
                        Try
                            mDatasetNum = Integer.Parse(strValue)
                        Catch ex As System.Exception
                            mDatasetNum = 0
                        End Try
                    End If

                    If .RetrieveValueForParameter("P", strValue) Then mParameterFilePath = strValue

                    If .RetrieveValueForParameter("Q", strValue) Then mQuietMode = True
                End With

                Return True
            End If

        Catch ex As System.Exception
            If mQuietMode Then
                Throw New System.Exception("Error parsing the command line parameters", ex)
            Else
                MsgBox("Error parsing the command line parameters: " & ControlChars.NewLine & ex.Message, MsgBoxStyle.Exclamation Or MsgBoxStyle.OKOnly, "Error")
            End If
        End Try

    End Function

    Private Sub ShowProgramHelp()

        Dim strSyntax As String
		Try

			strSyntax = "This program will read a text file with mass and intensity data for MS spectra and determine the LC-MS features present using UMCCreation.dll." & ControlChars.NewLine & ControlChars.NewLine
			strSyntax &= "Program syntax:" & ControlChars.NewLine & System.IO.Path.GetFileName(System.Reflection.Assembly.GetExecutingAssembly().Location)
			strSyntax &= " /I:InputFileName.txt [/O:OutputFileName.txt]" & ControlChars.NewLine & ControlChars.NewLine

			strSyntax &= "The input file name is required. If the filename contains spaces, then surround it with double quotes. " & _
						 "If the output file name is not supplied, then it will be auto-generated as InputFileName_Features.txt. " & _
						 "This program will look for a .Ini file with the same name as the input file, but with extension .Ini  " & _
						 "If found, then it reads the .Ini file to load the weight values to use when clustering the data." & ControlChars.NewLine & ControlChars.NewLine

			strSyntax &= "Program written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA) in 2007" & ControlChars.NewLine
			strSyntax &= "Copyright 2007, Battelle Memorial Institute.  All Rights Reserved." & ControlChars.NewLine & ControlChars.NewLine

			strSyntax &= "E-mail: matthew.monroe@pnl.gov or matt@alchemistmatt.com" & ControlChars.NewLine
			strSyntax &= "Website: http://ncrr.pnl.gov/ or http://www.sysbio.org/resources/staff/"

			Console.WriteLine(strSyntax)

		Catch ex As System.Exception
			Console.WriteLine("Error displaying the program syntax: " & ex.Message)
		End Try

    End Sub

End Module
