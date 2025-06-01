## Prerequisites ()
>=Python 3.8.10
https://learn.microsoft.com/en-US/cpp/windows/latest-supported-vc-redist?view=msvc-170#visual-studio-2015-2017-2019-and-2022   - for exe

### Usage
main.py - текущая версия программы (GUI-QT5) выбора любых данный из файлов (csv, gz, txt)
grath_matplot.py - модуль построения графиков на matplotlib



1.  RBMK_naladka.py(grath_from_tg_static_8.py - в предыдущих версиях) - текущая версия программы для работы на РБМК -> скачки, смещения и т.д.(GUI-tkinter) 
    
    grath_from_tg_static_9.py - аналог main.py только под(GUI-tkinter), для Xp - последняя версия питона 3.4.2, PyQT5 ставиться на 3.5 и выше
    grath_from_tg_static_9(Qt4).py - аналог read_csv_3_any только под(Qt4) через anaconda, для Xp - последняя версия питона 3.4.2, PyQT5 ставиться на 3.5 и выше


### help
2.  Создание exe (второй вариант)
    1. https://stackoverflow.com/questions/41570359/how-can-i-convert-a-py-to-exe-for-python
    1. setup.py - для создания exe через cx_Freeze
    2. в командной строке         python setup.py build


3. конвертирование ui в py (через терминал) ------ python -m PyQt5.uic.pyuic -x qt_des.ui -o qt_des.py


4. PyCharm 
Help -> Run->Edit Configurations -> Emulate terminal in output console* (для расшифровки ошибок  pyqt5)



TODO
1. Узнать сколько памяти использует программа (memory_profiler и psutil)
2. перейти на timeit
3. добавить логирование ()
4. Matplotlib заменить на PyQtGrath
5. Сделать аналогичную версию на С++
6. Сделать серверную службу
7. Сделать десктопную версию
8. Сделать мобильную версию
9. Добавить тесты
10. Поиск в списке сигналов (Qtablewidget) https://stackoverflow.com/questions/51613638/highlight-search-results-in-qtablewidgetselect-and-highlight-that-text-or-chara
11. Портировать на Linux(Astra и Ubuntu)
12. добавить The QProgressBar
13. 📦src
 ┣ 📂config
 ┃ ┣ 📜__init__.py
 ┃ ┗ 📜config.py
 ┣ 📂ui
 ┃ ┣ 📜__init__.py
 ┃ ┣ 📜grath_matplot.py
 ┃ ┣ 📜MainWindowUI.py
 ┣ 📂utils
 ┃ ┣ 📜__init__.py
 ┃ ┗ 📜helpers.py
 ┣ 📂logic
 ┃ ┣ 📜__init__.py
 ┃ ┣ 📜main_logic.py
 ┃ ┣ 📜data_processing.py
 ┃ ┗ 📜network_operations.py
 ┣ 📂logs
 ┃ ┗ 📜app.log
 ┣ 📂tests
 ┃ ┣ 📜__init__.py
 ┃ ┣ 📜test_main.py
 ┃ ┗ 📜test_logic.py
 ┣ 📂scripts
 ┃ ┣ 📜build.py
 ┃ ┣ 📜deploy.py
 ┣ 📂assets
 ┃ ┣ 📂images
 ┃ ┣ 📂fonts
 ┣ 📂docs
 ┃ ┣ 📜index.md
 ┃ ┣ 📜installation.md
 ┣ 📜main.py
 ┣ 📜MainWindowLogic.py
 ┣ 📜logger.py
 ┣ 📜README.md
 ┣ 📜requirements.txt
 ┗ 📜.gitignore

help:
python -m pip  install -r requirements.txt --no-index -f d:\\temp\\python_Library

### Ugage (Windows)
1. source build.bat -> создание текущеко билда