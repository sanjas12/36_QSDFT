# plot_test_data.py
import pandas as pd
import matplotlib.pyplot as plt

# Чтение и обработка данных (как в предыдущем коде)
with open('data_in.csv', 'r') as f:
    lines = f.readlines()

table1_start = None
table2_start = None
for i, line in enumerate(lines):
    if "Table 1" in line:
        table1_start = i + 1
    elif "Table 2" in line:
        table2_start = i + 1

# Чтение данных в DataFrame
df1 = pd.read_csv('data_in.csv', 
                 skiprows=table1_start, 
                 nrows=table2_start - table1_start - 2,
                 sep='\t',
                 engine='python')

df2 = pd.read_csv('data_in.csv', 
                 skiprows=table2_start, 
                 sep='\t',
                 engine='python')

# Очистка данных
df1 = df1.dropna(axis=1, how='all')
df2 = df2.dropna(axis=1, how='all')

# Присвоение имен столбцам
columns = ['Pore_size', 'Pressure', 'Pore_size_2', 'Pressure_2', 'Pore_size_3', 'Pressure_3']
df1.columns = columns
df2.columns = columns

# Создаем фигуру с несколькими графиками
plt.figure(figsize=(15, 10))

# График 1: Данные для азота (Table 1)
plt.subplot(2, 2, 1)
plt.plot(df1['Pore_size'], df1['Pressure'], 'b-', label='Set 1')
plt.plot(df1['Pore_size_2'], df1['Pressure_2'], 'g-', label='Set 2')
plt.plot(df1['Pore_size_3'], df1['Pressure_3'], 'r-', label='Set 3')
plt.title('Nitrogen Adsorption at 77.4K')
plt.xlabel('Pore size (Å)')
plt.ylabel('Filling pressure (P/P0)')
plt.yscale('log')
plt.legend()
plt.grid(True)

# График 2: Данные для аргона (Table 2)
plt.subplot(2, 2, 2)
plt.plot(df2['Pore_size'], df2['Pressure'], 'b-', label='Set 1')
plt.plot(df2['Pore_size_2'], df2['Pressure_2'], 'g-', label='Set 2')
plt.plot(df2['Pore_size_3'], df2['Pressure_3'], 'r-', label='Set 3')
plt.title('Argon Adsorption at 87.3K')
plt.xlabel('Pore size (Å)')
plt.ylabel('Filling pressure (P/P0)')
plt.yscale('log')
plt.legend()
plt.grid(True)

# График 3: Сравнение первых наборов данных
plt.subplot(2, 2, 3)
plt.plot(df1['Pore_size'], df1['Pressure'], 'b-', label='N2 at 77.4K')
plt.plot(df2['Pore_size'], df2['Pressure'], 'r--', label='Ar at 87.3K')
plt.title('Comparison of First Data Sets')
plt.xlabel('Pore size (Å)')
plt.ylabel('Filling pressure (P/P0)')
plt.yscale('log')
plt.legend()
plt.grid(True)

# График 4: Сравнение всех данных (усреднённое)
plt.subplot(2, 2, 4)
# Объединяем все данные по порам для каждого газа
all_n2 = pd.concat([df1['Pore_size'], df1['Pore_size_2'], df1['Pore_size_3']])
all_n2_pressure = pd.concat([df1['Pressure'], df1['Pressure_2'], df1['Pressure_3']])
all_ar = pd.concat([df2['Pore_size'], df2['Pore_size_2'], df2['Pore_size_3']])
all_ar_pressure = pd.concat([df2['Pressure'], df2['Pressure_2'], df2['Pressure_3']])

plt.scatter(all_n2, all_n2_pressure, c='blue', s=10, alpha=0.5, label='N2 at 77.4K')
plt.scatter(all_ar, all_ar_pressure, c='red', s=10, alpha=0.5, label='Ar at 87.3K')
plt.title('All Data Comparison')
plt.xlabel('Pore size (Å)')
plt.ylabel('Filling pressure (P/P0)')
plt.yscale('log')
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.show()