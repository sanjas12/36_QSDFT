# plot_test_data_class.py
import pandas as pd
import matplotlib.pyplot as plt

class AdsorptionDataPlotter:
    def __init__(self, file_path='data_in.csv'):
        """
        Initialize the plotter with data file path.
        
        Args:
            file_path (str): Path to the input CSV file containing adsorption data
        """
        self.file_path = file_path
        self.df1 = None  # Table 1 data (Nitrogen)
        self.df2 = None  # Table 2 data (Argon)
        
    def load_data(self):
        """
        Load and process data from the input file.
        Splits data into two tables and cleans them.
        """
        with open(self.file_path, 'r') as f:
            lines = f.readlines()

        # Find table boundaries
        table1_start = None
        table2_start = None
        for i, line in enumerate(lines):
            if "Table 1" in line:
                table1_start = i + 1
            elif "Table 2" in line:
                table2_start = i + 1

        # Read data into DataFrames
        self.df1 = pd.read_csv(self.file_path, 
                             skiprows=table1_start, 
                             nrows=table2_start - table1_start - 2,
                             sep='\t',
                             engine='python')

        self.df2 = pd.read_csv(self.file_path, 
                             skiprows=table2_start, 
                             sep='\t',
                             engine='python')

        # Clean data
        self.df1 = self.df1.dropna(axis=1, how='all')
        self.df2 = self.df2.dropna(axis=1, how='all')

        # Set column names
        columns = ['Pore_size', 'Pressure', 'Pore_size_2', 'Pressure_2', 'Pore_size_3', 'Pressure_3']
        self.df1.columns = columns
        self.df2.columns = columns
        
    def plot_all(self, figsize=(15, 10)):
        """
        Create a figure with 4 subplots showing different views of the data.
        
        Args:
            figsize (tuple): Size of the output figure
        """
        if self.df1 is None or self.df2 is None:
            self.load_data()
            
        plt.figure(figsize=figsize)
        
        # Plot 1: Nitrogen data
        self._plot_single_gas(1, self.df1, 'Nitrogen Adsorption at 77.4K')
        
        # Plot 2: Argon data
        self._plot_single_gas(2, self.df2, 'Argon Adsorption at 87.3K')
        
        # Plot 3: First set comparison
        self._plot_comparison(3, 
                            (self.df1['Pore_size'], self.df1['Pressure'], 'N2 at 77.4K', 'b-'),
                            (self.df2['Pore_size'], self.df2['Pressure'], 'Ar at 87.3K', 'r--'),
                            'Comparison of First Data Sets')
        
        # Plot 4: All data comparison
        self._plot_all_data_comparison(4)
        
        plt.tight_layout()
        plt.show()
        
    def _plot_single_gas(self, subplot_num, df, title):
        """
        Helper method to plot data for a single gas.
        """
        plt.subplot(2, 2, subplot_num)
        plt.plot(df['Pore_size'], df['Pressure'], 'b-', label='Set 1')
        plt.plot(df['Pore_size_2'], df['Pressure_2'], 'g-', label='Set 2')
        plt.plot(df['Pore_size_3'], df['Pressure_3'], 'r-', label='Set 3')
        plt.title(title)
        plt.xlabel('Pore size (Å)')
        plt.ylabel('Filling pressure (P/P0)')
        plt.yscale('log')
        plt.legend()
        plt.grid(True)
        
    def _plot_comparison(self, subplot_num, data1, data2, title):
        """
        Helper method to plot comparison between two datasets.
        """
        plt.subplot(2, 2, subplot_num)
        x1, y1, label1, style1 = data1
        x2, y2, label2, style2 = data2
        plt.plot(x1, y1, style1, label=label1)
        plt.plot(x2, y2, style2, label=label2)
        plt.title(title)
        plt.xlabel('Pore size (Å)')
        plt.ylabel('Filling pressure (P/P0)')
        plt.yscale('log')
        plt.legend()
        plt.grid(True)
        
    def _plot_all_data_comparison(self, subplot_num):
        """
        Helper method to plot all data points comparison.
        """
        # Combine all data for each gas
        all_n2 = pd.concat([self.df1['Pore_size'], self.df1['Pore_size_2'], self.df1['Pore_size_3']])
        all_n2_pressure = pd.concat([self.df1['Pressure'], self.df1['Pressure_2'], self.df1['Pressure_3']])
        all_ar = pd.concat([self.df2['Pore_size'], self.df2['Pore_size_2'], self.df2['Pore_size_3']])
        all_ar_pressure = pd.concat([self.df2['Pressure'], self.df2['Pressure_2'], self.df2['Pressure_3']])

        plt.subplot(2, 2, subplot_num)
        plt.scatter(all_n2, all_n2_pressure, c='blue', s=10, alpha=0.5, label='N2 at 77.4K')
        plt.scatter(all_ar, all_ar_pressure, c='red', s=10, alpha=0.5, label='Ar at 87.3K')
        plt.title('All Data Comparison')
        plt.xlabel('Pore size (Å)')
        plt.ylabel('Filling pressure (P/P0)')
        plt.yscale('log')
        plt.legend()
        plt.grid(True)


# Example usage:
if __name__ == "__main__":
    plotter = AdsorptionDataPlotter('data_in.csv')
    plotter.plot_all()