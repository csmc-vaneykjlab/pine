set "PATH=%PATH%;C:\Windows\System32\downlevel;"
env\Scripts\activate.ps1
pyinstaller --onedir --clean --hidden-import igraph.vendor.texttable .\pine\pine.py
deactivate