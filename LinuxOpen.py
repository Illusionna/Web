import os
import webbrowser

os.system('clear')
os.system('make html')

webbrowser.open(
    url = os.getcwd() + os.sep + './build/html/index.html',
    new = 0,
    autoraise = True
)