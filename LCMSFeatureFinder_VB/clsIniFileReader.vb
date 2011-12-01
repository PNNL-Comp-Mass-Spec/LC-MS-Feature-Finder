Option Strict On

' This class can read a .Ini file containing configuration settings
' One can obtain the value for a given setting in a given section using GetSetting()
' One can also update values or add new values and sections
' To save an updated file, call .WriteIniFile() or .WriteIniFile(strIniFilePath)
' Note that Section names and Key names are not case sensitive (though the original case is preserved when writing out the data to a file)
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

Public Class clsIniFileReader

    Public Sub New(ByVal strIniFilePath As String)
        ReadIniFile(strIniFilePath)
    End Sub

#Region "Constants and Enums"
    Protected Const COMMENT_LINE_START_CHAR As Char = ";"c
    Protected Const SECTION_LINE_START_CHAR As Char = "["c
    Protected Const SECTION_LINE_END_CHAR As Char = "]"c

    Protected Const INITIAL_SECTION_COUNT As Integer = 4
    Protected Const INITIAL_KEY_COUNT_PER_SECTION As Integer = 2
#End Region

#Region "Structures"

    Protected Structure udtKeyValueListType
        Public Count As Integer
        Public NamesLowerCase() As String
        Public NamesOrignalCase() As String
        Public Values() As String
    End Structure

#End Region

#Region "Classwide Variables"
    Protected mIniFilePath As String
    Protected mDataLoaded As Boolean

    Protected mSectionCount As Integer
    Protected mSectionsLowerCase() As String
    Protected mSectionsOriginalCase() As String
    Protected mKeysBySection() As udtKeyValueListType
#End Region

#Region "Properties"
    Public ReadOnly Property DataLoaded() As Boolean
        Get
            Return mDataLoaded
        End Get
    End Property

    Public ReadOnly Property IniFilePath() As String
        Get
            Return mIniFilePath
        End Get
    End Property

    Public ReadOnly Property SectionNames() As String()
        Get
            Return mSectionsOriginalCase
        End Get
    End Property
#End Region

    Protected Sub ClearSectionData()
        Dim intIndex As Integer

        mSectionCount = 0

        ReDim mSectionsLowerCase(INITIAL_SECTION_COUNT - 1)
        ReDim mSectionsOriginalCase(INITIAL_SECTION_COUNT - 1)
        ReDim mKeysBySection(INITIAL_SECTION_COUNT - 1)

        For intIndex = 0 To mKeysBySection.Length - 1
            With mKeysBySection(intIndex)
                .Count = 0
                ReDim .NamesLowerCase(INITIAL_KEY_COUNT_PER_SECTION - 1)
                ReDim .NamesOrignalCase(INITIAL_KEY_COUNT_PER_SECTION - 1)
                ReDim .Values(INITIAL_KEY_COUNT_PER_SECTION - 1)
            End With
        Next

        mDataLoaded = False
    End Sub

    Protected Sub ExpandSectionDataIfNeeded(ByVal intLengthRequired As Integer)
        Dim intIndex As Integer
        Dim intStartIndex As Integer
        Dim intNewLength As Integer

        If intLengthRequired > mSectionsLowerCase.Length Then
            intStartIndex = mSectionsLowerCase.Length

            intNewLength = mSectionsLowerCase.Length * 2
            Do While intNewLength < intLengthRequired
                intNewLength *= 2
            Loop

            ReDim Preserve mSectionsLowerCase(intNewLength - 1)
            ReDim Preserve mSectionsOriginalCase(intNewLength - 1)
            ReDim Preserve mKeysBySection(intNewLength - 1)

            For intIndex = intStartIndex To mKeysBySection.Length - 1
                With mKeysBySection(intIndex)
                    .Count = 0
                    ReDim .NamesLowerCase(INITIAL_KEY_COUNT_PER_SECTION - 1)
                    ReDim .NamesOrignalCase(INITIAL_KEY_COUNT_PER_SECTION - 1)
                    ReDim .Values(INITIAL_KEY_COUNT_PER_SECTION - 1)
                End With
            Next
        End If
    End Sub

    Protected Sub ExpandKeyDataIfNeeded(ByVal intSectionIndex As Integer, ByVal intLengthRequired As Integer)
		Dim intNewLength As Integer

        If intSectionIndex >= 0 And intSectionIndex < mSectionsLowerCase.Length Then
            If intLengthRequired > mKeysBySection(intSectionIndex).NamesLowerCase.Length Then
                With mKeysBySection(intSectionIndex)
                    intNewLength = .NamesLowerCase.Length * 2
                    Do While intNewLength < intLengthRequired
                        intNewLength *= 2
                    Loop
                    ReDim Preserve .NamesLowerCase(intNewLength - 1)
                    ReDim Preserve .NamesOrignalCase(intNewLength - 1)
                    ReDim Preserve .Values(intNewLength - 1)
                End With
            End If
        End If

    End Sub

    Protected Function GetSectionIndex(ByVal strSectionName As String, ByVal blnAddIfMissing As Boolean) As Integer
        ' Looks for strSectionName in mSectionsLowerCase
        ' If not found, then adds it


        Dim intSectionIndex As Integer

        If mSectionsLowerCase Is Nothing Then
            ClearSectionData()
        End If

        intSectionIndex = GetSectionIndex(strSectionName)
        If intSectionIndex < 0 AndAlso blnAddIfMissing Then
            ExpandSectionDataIfNeeded(mSectionCount + 1)

            intSectionIndex = mSectionCount
            mSectionCount += 1

            mSectionsLowerCase(intSectionIndex) = strSectionName.ToLower
            mSectionsOriginalCase(intSectionIndex) = String.Copy(strSectionName)

            mKeysBySection(intSectionIndex).Count = 0
        End If

        Return intSectionIndex

    End Function

    Protected Function GetSettingConfirmMatch(ByVal strSectionName As String, ByVal strKeyName As String, ByVal strDefaultValue As String, ByRef blnMatchFound As Boolean) As String
        Dim intSectionIndex As Integer
        Dim intIndexMatch As Integer

        Dim strValue As String

        strValue = String.Copy(strDefaultValue)
        blnMatchFound = False

        intSectionIndex = GetSectionIndex(strSectionName)
        If intSectionIndex >= 0 Then
            With mKeysBySection(intSectionIndex)
                intIndexMatch = Array.IndexOf(.NamesLowerCase, strKeyName.ToLower, 0, .Count)

                If intIndexMatch >= 0 Then
                    strValue = .Values(intIndexMatch)
                    blnMatchFound = True
                End If
            End With
        End If

        Return strValue
    End Function

    Public Function GetSetting(ByVal strSectionName As String, ByVal strKeyName As String, ByVal strDefaultValue As String) As String
        Return GetSettingConfirmMatch(strSectionName, strKeyName, strDefaultValue, False)
    End Function

    Public Function GetSetting(ByVal strSectionName As String, ByVal strKeyName As String, ByVal intDefaultValue As Integer) As Integer
        Dim blnMatchFound As Boolean
        Dim strValue As String
        Dim intValue As Integer

        intValue = intDefaultValue
        strValue = GetSettingConfirmMatch(strSectionName, strKeyName, intDefaultValue.ToString, blnMatchFound)

        If blnMatchFound Then
            Try
                If IsNumber(strValue) Then
                    intValue = CInt(strValue)
                End If
            Catch ex As System.Exception
                ' Ignore errors here
            End Try
        End If

        Return intValue

    End Function

    Public Function GetSetting(ByVal strSectionName As String, ByVal strKeyName As String, ByVal sngDefaultValue As Single) As Single
        Dim blnMatchFound As Boolean
        Dim strValue As String
        Dim sngValue As Single

        sngValue = sngDefaultValue
        strValue = GetSettingConfirmMatch(strSectionName, strKeyName, sngDefaultValue.ToString, blnMatchFound)

        If blnMatchFound Then
            Try
                If IsNumber(strValue) Then
                    sngValue = CSng(strValue)
                End If
            Catch ex As System.Exception
                ' Ignore errors here
            End Try
        End If

        Return sngValue

    End Function

    Public Function GetSetting(ByVal strSectionName As String, ByVal strKeyName As String, ByVal dblDefaultValue As Double) As Double
        Dim blnMatchFound As Boolean
        Dim strValue As String
        Dim dblValue As Double

        dblValue = dblDefaultValue
        strValue = GetSettingConfirmMatch(strSectionName, strKeyName, dblDefaultValue.ToString, blnMatchFound)

        If blnMatchFound Then
            Try
                If IsNumber(strValue) Then
                    dblValue = CDbl(strValue)
                End If
            Catch ex As System.Exception
                ' Ignore errors here
            End Try
        End If

        Return dblValue

    End Function

    Public Function GetSetting(ByVal strSectionName As String, ByVal strKeyName As String, ByVal blnDefaultValue As Boolean) As Boolean
        Dim blnMatchFound As Boolean
        Dim strValue As String
        Dim blnValue As Boolean

        blnValue = blnDefaultValue
        strValue = GetSettingConfirmMatch(strSectionName, strKeyName, blnDefaultValue.ToString, blnMatchFound)

        If blnMatchFound Then
            Try
                blnValue = CBool(strValue)
            Catch ex As System.Exception
                ' Ignore errors here
            End Try
        End If

        Return blnValue

    End Function

    Protected Function IsNumber(ByVal strValue As String) As Boolean
        Dim objFormatProvider As System.Globalization.NumberFormatInfo
        Try
            Return Double.TryParse(strValue, Globalization.NumberStyles.Any, objFormatProvider, 0)
        Catch ex As System.Exception
            Return False
        End Try
    End Function

    Public Function ReadIniFile(ByVal strIniFilePath As String) As Boolean
        ' Reads the settings from the .Ini file specifed by strIniFilePath
        ' Stores the values in mSectionsLowerCase, mSectionsOriginalCase and mKeysBySection

        Dim srInFile As System.IO.StreamReader

        Dim strLineIn As String

        Dim strCurrentSection As String
        Dim intCurrentSectionIndex As Integer

        Dim strKeyName As String
        Dim strValue As String

        Dim intEqualsLoc As Integer
        Dim blnSuccess As Boolean

        Try
            If Not System.IO.File.Exists(strIniFilePath) Then
                Return False
            End If

            ' Assume success for now
            blnSuccess = True

            ClearSectionData()
            mIniFilePath = String.Empty

            srInFile = New System.IO.StreamReader(New System.IO.FileStream(strIniFilePath, IO.FileMode.Open, IO.FileAccess.Read, IO.FileShare.Read))

            mIniFilePath = String.Copy(strIniFilePath)
            strCurrentSection = String.Empty

            Do While srInFile.Peek >= 0
                strLineIn = srInFile.ReadLine

                If Not strLineIn Is Nothing AndAlso strLineIn.Length > 0 Then
                    strLineIn = strLineIn.Trim

                    If strLineIn.StartsWith(COMMENT_LINE_START_CHAR) Then
                        ' Comment line; skip it
                    Else
                        If strLineIn.StartsWith(SECTION_LINE_START_CHAR) AndAlso strLineIn.EndsWith(SECTION_LINE_END_CHAR) Then
                            ' Section name
                            strCurrentSection = strLineIn.Substring(1, strLineIn.Length - 2)
                            intCurrentSectionIndex = GetSectionIndex(strCurrentSection, True)
                        Else
                            ' Look for an Equals sign
                            intEqualsLoc = strLineIn.IndexOf("="c)

                            If intEqualsLoc > 0 Then
                                strKeyName = strLineIn.Substring(0, intEqualsLoc)
                                strValue = strLineIn.Substring(intEqualsLoc + 1)

                                StoreKeyAndValue(intCurrentSectionIndex, strKeyName, strValue)
                            End If
                        End If
                    End If
                End If
            Loop

            If mSectionCount > 0 Then
                mDataLoaded = True
            End If

            If mSectionCount < mSectionsLowerCase.Length Then
                ReDim Preserve mSectionsLowerCase(mSectionCount - 1)
                ReDim Preserve mSectionsOriginalCase(mSectionCount - 1)
                ReDim Preserve mKeysBySection(mSectionCount - 1)
            End If
        Catch ex As System.Exception
            Console.WriteLine("Error in clsIniFileReader->ReadIniFile: " & ex.Message)
            blnSuccess = False
        Finally
            If Not srInFile Is Nothing Then
                srInFile.Close()
            End If
        End Try

        Return blnSuccess

    End Function

    Protected Function GetSectionIndex(ByVal strSectionName As String) As Integer
        Dim intSectionIndex As Integer

        If mSectionsLowerCase Is Nothing Then
            intSectionIndex = -1
        Else
            intSectionIndex = Array.IndexOf(mSectionsLowerCase, strSectionName.ToLower)
        End If

        Return intSectionIndex
    End Function

    Public Sub SaveSetting(ByVal strSectionName As String, ByVal strKeyName As String, ByVal strValue As String)
        StoreKeyAndValue(strSectionName, strKeyName, strValue)
    End Sub

    Public Sub SaveSetting(ByVal strSectionName As String, ByVal strKeyName As String, ByVal intValue As Integer)
        StoreKeyAndValue(strSectionName, strKeyName, intValue.ToString)
    End Sub

    Public Sub SaveSetting(ByVal strSectionName As String, ByVal strKeyName As String, ByVal dblValue As Double)
        StoreKeyAndValue(strSectionName, strKeyName, dblValue.ToString)
    End Sub

    Public Sub SaveSetting(ByVal strSectionName As String, ByVal strKeyName As String, ByVal sngValue As Single)
        StoreKeyAndValue(strSectionName, strKeyName, sngValue.ToString)
    End Sub

    Public Sub SaveSetting(ByVal strSectionName As String, ByVal strKeyName As String, ByVal blnValue As Boolean)
        StoreKeyAndValue(strSectionName, strKeyName, blnValue.ToString)
    End Sub

    Protected Sub StoreKeyAndValue(ByVal strSectionName As String, ByVal strKeyName As String, ByVal strValue As String)
        Dim intSectionIndex As Integer

        intSectionIndex = GetSectionIndex(strSectionName, True)
        If intSectionIndex >= 0 Then
            StoreKeyAndValue(intSectionIndex, strKeyName, strValue)
        End If
    End Sub

    Protected Sub StoreKeyAndValue(ByVal intSectionIndex As Integer, ByVal strKeyName As String, ByVal strValue As String)
        Dim intIndexMatch As Integer

        If intSectionIndex >= 0 And intSectionIndex < mKeysBySection.Length Then
            With mKeysBySection(intSectionIndex)
                intIndexMatch = Array.IndexOf(.NamesLowerCase, strKeyName.ToLower, 0, .Count)

                If intIndexMatch >= 0 Then
                    .Values(intIndexMatch) = strValue
                Else
                    ExpandKeyDataIfNeeded(intSectionIndex, .Count + 1)
                    .NamesOrignalCase(.Count) = String.Copy(strKeyName)
                    .NamesLowerCase(.Count) = strKeyName.ToLower
                    .Values(.Count) = strValue
                    .Count += 1
                End If
            End With
        End If

    End Sub

    Public Function WriteIniFile() As Boolean
        If mDataLoaded Then
            Return WriteIniFile(mIniFilePath)
        Else
            Return False
        End If
    End Function

    Public Function WriteIniFile(ByVal strIniFilePath As String) As Boolean
        Dim srOutFile As System.IO.StreamWriter

        Dim blnSuccess As Boolean

        Dim intSectionIndex As Integer
        Dim intKeyIndex As Integer
        Dim strValue As String

        If mDataLoaded Then
            Try
                ' Assume success for now
                blnSuccess = True

                srOutFile = New System.IO.StreamWriter(New System.IO.FileStream(strIniFilePath, IO.FileMode.Create, IO.FileAccess.Write, IO.FileShare.Read))

                ' Update mIniFilePath
                mIniFilePath = String.Copy(strIniFilePath)

                For intSectionIndex = 0 To mSectionCount - 1
                    srOutFile.WriteLine(SECTION_LINE_START_CHAR & mSectionsOriginalCase(intSectionIndex) & SECTION_LINE_END_CHAR)

                    With mKeysBySection(intSectionIndex)
                        For intKeyIndex = 0 To .Count - 1
                            If Not .NamesOrignalCase(intKeyIndex) Is Nothing Then
                                strValue = String.Copy(.Values(intKeyIndex))
                                If strValue Is Nothing Then strValue = String.Empty

                                srOutFile.WriteLine(.NamesOrignalCase(intKeyIndex) & "=" & strValue)
                            End If
                        Next intKeyIndex
                    End With

                Next intSectionIndex

            Catch ex As System.Exception
                blnSuccess = False
            Finally
                If Not srOutFile Is Nothing Then
                    srOutFile.Close()
                End If
            End Try

        Else
            blnSuccess = False
        End If

        Return blnSuccess
    End Function

End Class

