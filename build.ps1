set "PATH=%PATH%;C:\Windows\System32\downlevel;"
env\Scripts\activate.ps1
pyinstaller --onedir --clean --noconfirm --hidden-import igraph.vendor.texttable --exclude-module matplotlib .\pine\pine.py
deactivate