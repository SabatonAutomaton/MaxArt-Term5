#pragma comment(linker, \
    "\"/manifestdependency:type='win32' \
    name='Microsoft.Windows.Common-Controls' version='6.0.0.0' \
    processorArchitecture='*' publicKeyToken='6595b64144ccf1df' language='*'\"")

#include <windows.h>

HWND hWndTask1, hWndExit;
HANDLE hThread;
wchar_t Info[256] = L"Задание 1 - текущий час:\n(загрузка)";
int CurrentHour( )
{
   SYSTEMTIME systemTime;
   SYSTEMTIME localSystemTime;
   GetSystemTime( &systemTime );
   SystemTimeToTzSpecificLocalTime( NULL, &systemTime, &localSystemTime );
   return localSystemTime.wHour;
}

DWORD WINAPI ThreadFunc( LPVOID lpParam )
{
   typedef int( *ImportFunciton )(WCHAR*);
   ImportFunciton DLLInfo;
   HINSTANCE hinstLib = LoadLibrary(TEXT("CurrentHour.dll" ) );
   DLLInfo = ( ImportFunciton )GetProcAddress(hinstLib,"CurrentHour" );
   
   int hour = DLLInfo( Info);
   //wchar_t buf[256];
   //wsprintfW( buf, L"Задание 1 - текущий час: %d", hour );
   //wcscpy_s( Info, 256, buf );

   return 0;
}

LRESULT CALLBACK WindowFunc( HWND hWnd, UINT msg, WPARAM wParam, LPARAM lParam )
{
   PAINTSTRUCT ps;
   HDC hDC;

   switch ( msg )
   {
      case WM_CREATE:
         //hWndTask1 = CreateWindow( L"STATIC", L"Задание 1 - текущий час:\n(загрузка)",
         //   WS_VISIBLE | WS_CHILD,
         //   10, 10, 440, 40,
         //   hWnd, ( HMENU )1, NULL, NULL );
         //SendMessage( hWndTask1, WM_SETFONT, NULL, TRUE );

         hWndExit = CreateWindow( L"BUTTON", L"Выйти",
            WS_TABSTOP | WS_VISIBLE | WS_CHILD | BS_DEFPUSHBUTTON,
            18, 110, 250, 30,
            hWnd, ( HMENU )3, ( HINSTANCE )GetWindowLongPtr( hWnd, GWLP_HINSTANCE ), NULL );
         SendMessage( hWndExit, WM_SETFONT, NULL, TRUE );

         hThread = CreateThread( NULL, 0, ThreadFunc, hWnd, 0, NULL );
         if ( hThread == NULL ) {
            MessageBox( NULL, L"Невозможно запустить поток!", L"Ошибка!", MB_ICONEXCLAMATION | MB_OK );
            break;
         }
         CloseHandle( hThread );
         break;
      case WM_COMMAND: // Взаимодействие с элементами управления
         // Переход к обработчику выхода только если нажата кнопка "Выйти"
         if ( LOWORD( wParam ) != 3 )
            break;
      case WM_DESTROY:
         PostQuitMessage( 0 );
         break;
      case WM_PAINT:
         hDC = BeginPaint( hWnd, &ps );
         TextOut( hDC, 10, 10, Info, wcslen( Info ) );
         EndPaint( hWnd, &ps );
         break;
      default:
         return DefWindowProc( hWnd, msg, wParam, lParam );
   }
   return 0;

}

int WINAPI WinMain( HINSTANCE hThisInst, HINSTANCE hPrevInst, LPSTR str, int nWinMode )
{
   SetConsoleCP( 65001 );
   SetConsoleOutputCP( 65001 );
   MSG msg;
   WNDCLASS wcl;
   HWND hWnd;
   wcl.hInstance = hThisInst;
   wcl.lpszClassName = L"main";
   wcl.lpfnWndProc = WindowFunc;
   wcl.style = CS_HREDRAW;
   wcl.hIcon = LoadIcon( NULL, IDI_APPLICATION );
   wcl.hCursor = LoadCursor( NULL, IDC_ARROW );
   wcl.lpszMenuName = NULL;
   wcl.cbWndExtra = 0;
   wcl.hbrBackground = ( HBRUSH )GetStockObject( WHITE_BRUSH );
   RegisterClass( &wcl );

   hWnd = CreateWindow( L"main", L"РГЗ ПМ-34 Рубцов А.В.", WS_OVERLAPPEDWINDOW, 100, 100, 300, 200, HWND_DESKTOP, NULL, hThisInst, NULL );
   ShowWindow( hWnd, nWinMode );

   while ( GetMessage( &msg, NULL, 0, 0 ) )
   {
      TranslateMessage( &msg );
      DispatchMessage( &msg );

   }
   return msg.wParam;
}