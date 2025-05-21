# plot_test_data.py
import pandas as pd

# Чтение всего файла
with open('data_in.csv', 'r') as f:
    lines = f.readlines()

# Находим индексы начала таблиц
table1_start = None
table2_start = None
for i, line in enumerate(lines):
    if "Table 1" in line:
        table1_start = i + 1  # Следующая строка после заголовка
    elif "Table 2" in line:
        table2_start = i + 1  # Следующая строка после заголовка

# Читаем первую таблицу
df1 = pd.read_csv('data_in.csv', 
                 skiprows=table1_start,  # type: ignore
                 nrows=table2_start - table1_start - 2,  # -2 чтобы исключить пустую строку и заголовок Table 2 # type: ignore
                 sep='\t',
                 engine='python')

# Читаем вторую таблицу
df2 = pd.read_csv('data_in.csv', 
                 skiprows=table2_start,  # type: ignore
                 sep='\t',
                 engine='python')

# Удаляем возможные NaN столбцы (если есть пустые столбцы из-за табуляции)
df1 = df1.dropna(axis=1, how='all')
df2 = df2.dropna(axis=1, how='all')

# Переименовываем столбцы (если нужно)
df1.columns = ['Pore size(Å)', 'Filling pressre(P/P0)', 'Pore size(Å)_2', 'Filling pressre(P/P0)_2', 'Pore size(Å)_3', 'Filling pressre(P/P0)_3']
df2.columns = ['Pore size(Å)', 'Filling pressre(P/P0)', 'Pore size(Å)_2', 'Filling pressre(P/P0)_2', 'Pore size(Å)_3', 'Filling pressre(P/P0)_3']

print("DataFrame 1:")
print(df1.head())
print("\nDataFrame 2:")
print(df2.head())