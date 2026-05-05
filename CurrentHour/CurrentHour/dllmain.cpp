// dllmain.cpp : Определяет точку входа для приложения DLL.
#include "pch.h"
#include <windows.h>

extern "C" _declspec( dllexport ) int CurrentHour( WCHAR *InfoString )
{
   SYSTEMTIME systemTime;
   SYSTEMTIME localSystemTime;
   GetSystemTime( &systemTime );
   SystemTimeToTzSpecificLocalTime( NULL, &systemTime, &localSystemTime );
   wchar_t buf[256];
   wsprintfW( buf, L"Задание 1 - текущий час: %d", localSystemTime.wHour );
   wcscpy_s( InfoString, 256, buf );
   return 0;
}

