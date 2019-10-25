.\node_modules\.bin\electron-packager --overwrite . pine
Copy-Item "..\dist\changes_to_pine_final.exe" -Destination ".\pine-win32-x64\resources\pine_2.exe"
#robocopy ..\dist\pine_2 .\pine-win32-x64\resources\pine_2 /E