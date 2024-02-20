import os
import webbrowser

os.system('cls')
os.system('make html')

webbrowser.register(
    name = 'chrome',
    klass = None,
    instance = webbrowser.BackgroundBrowser(
        name = './Chromium/chrome.exe'
    )
)
browser = webbrowser.get('chrome')
browser.open(
    url = os.getcwd() + os.sep + './build/html/index.html',
    new = 0,
    autoraise = True
)
