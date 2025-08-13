"""
GermiTrack: Smart Analysis for Germination Data - v.1.4 DESKTOP

Professional desktop application for germination analysis

Authors: Pontes et al. (2025)
Developed by GOF/UFMS
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import math
import warnings
import os
import sys
from pathlib import Path
import tkinter as tk
from tkinter import filedialog, messagebox, ttk, scrolledtext

warnings.filterwarnings('ignore')

# ==========================================
# 1. Analysis and Export Functions (Modular)
# ==========================================

def save_germination_curves_data(data, treatment_name="Treatment", save_path=""):
    rep_columns = [col for col in data.columns if col.startswith('R')]
    if not rep_columns:
        rep_columns = data.columns[1:]

    if 'Day' in data.columns:
        days = data['Day']
    elif 'ti' in data.columns:
        days = data['ti']
    else:
        days = data.iloc[:, 0]

    cumulative_data = {}
    daily_data = {}

    for col in rep_columns[:6]:
        if col in data.columns:
            daily = data[col].fillna(0).values
            cumulative = np.cumsum(daily)
            cumulative_data[col] = cumulative
            daily_data[col] = daily

    cumulative_array = np.array(list(cumulative_data.values())).T if cumulative_data else np.array([])
    mean_cumulative = np.mean(cumulative_array, axis=1) if cumulative_array.size > 0 else np.array([])
    std_cumulative = np.std(cumulative_array, axis=1) if cumulative_array.size > 0 else np.array([])

    rep_data = data[rep_columns].fillna(0)
    mean_daily = rep_data.mean(axis=1)
    std_daily = rep_data.std(axis=1)

    df_save = pd.DataFrame({'Day': days})
    for col, cum in cumulative_data.items():
        df_save[f'Cumulative_{col}'] = cum
    if mean_cumulative.size > 0:
        df_save['Mean_Cumulative'] = mean_cumulative
        df_save['Std_Cumulative'] = std_cumulative
    for col, daily in daily_data.items():
        df_save[f'Daily_{col}'] = daily
    df_save['Mean_Daily'] = mean_daily
    df_save['Std_Daily'] = std_daily

    if save_path:
        excel_file = f'{save_path}/GermiTrack_Germination_Curves_{treatment_name}.xlsx'
        csv_file = f'{save_path}/GermiTrack_Germination_Curves_{treatment_name}.csv'
    else:
        excel_file = f'GermiTrack_Germination_Curves_{treatment_name}.xlsx'
        csv_file = f'GermiTrack_Germination_Curves_{treatment_name}.csv'

    df_save.to_excel(excel_file, index=False)
    df_save.to_csv(csv_file, index=False)

    return excel_file, csv_file

class GermiTrackAnalyzer:
    def __init__(self, data, total_seeds=25):
        self.data = data
        self.total_seeds = total_seeds

    def calculate_all_germination_parameters(self, treatment_data):
        results = []
        rep_columns = [col for col in treatment_data.columns if col.startswith('R') or 'Rep' in col]
        if not rep_columns:
            rep_columns = treatment_data.columns[1:]
        if 'Day' in treatment_data.columns:
            days = treatment_data['Day'].values
        elif 'ti' in treatment_data.columns:
            days = treatment_data['ti'].values
        else:
            days = treatment_data.iloc[:, 0].values
        days = np.array(days)
        for i, rep_col in enumerate(rep_columns):
            if rep_col in treatment_data.columns:
                rep_data = treatment_data[rep_col].fillna(0).values
                total_germinated = np.sum(rep_data)
                if total_germinated > 0:
                    germinability = (total_germinated / self.total_seeds) * 100
                    mgt = np.sum(days * rep_data) / total_germinated
                    mgr = 1 / mgt if mgt > 0 else 0
                    variance = np.sum(rep_data * (days - mgt) ** 2) / (total_germinated - 1) if total_germinated > 1 else 0
                    std_dev = np.sqrt(variance)
                    cv = (std_dev / mgt * 100) if mgt > 0 else 0
                    speed_maguire = np.sum(rep_data / days) if np.all(days > 0) else 0
                    frequencies = rep_data / total_germinated
                    frequencies = frequencies[frequencies > 0]
                    uncertainty = -np.sum(frequencies * np.log2(frequencies)) if len(frequencies) > 0 else 0
                    numerator = sum(self._combination(ni, 2) for ni in rep_data if ni >= 2)
                    denominator = self._combination(total_germinated, 2)
                    synchrony = numerator / denominator if denominator > 0 else 0
                    t50 = self._calculate_tx(days, rep_data, total_germinated, 0.50)
                else:
                    germinability = mgt = mgr = variance = std_dev = cv = 0
                    speed_maguire = uncertainty = synchrony = t50 = 0
                results.append({
                    'Replicate': f'R{i+1}',
                    'Total_Seeds': self.total_seeds,
                    'Germinated_Seeds': int(total_germinated),
                    'G_Germinability_%': round(germinability, 2),
                    'MT_Mean_Germination_Time': round(mgt, 2),
                    'CVt_Coefficient_Variation_%': round(cv, 2),
                    'MR_Mean_Germination_Rate': round(mgr, 4),
                    'U_Uncertainty_Index': round(uncertainty, 3),
                    'Z_Synchrony_Index': round(synchrony, 3),
                    'Variance_Germination_Time': round(variance, 2),
                    'Standard_Deviation': round(std_dev, 2),
                    'Speed_Maguire_Index': round(speed_maguire, 2),
                    'T50_Time_50%': round(t50, 2),
                    'ArcSin_Transformation': round(np.arcsin(np.sqrt(germinability / 100)) * (180 / np.pi), 2) if germinability > 0 else 0
                })
        return pd.DataFrame(results)

    def _combination(self, n, r):
        if n < r or n < 0 or r < 0:
            return 0
        try:
            return math.factorial(n) / (math.factorial(r) * math.factorial(n - r))
        except:
            return 0

    def _calculate_tx(self, days, rep_data, total_germinated, percentage):
        if total_germinated == 0:
            return 0
        cumulative = np.cumsum(rep_data)
        target = total_germinated * percentage
        for i, cum_germ in enumerate(cumulative):
            if cum_germ >= target:
                if i == 0:
                    return days[0]
                else:
                    x1, y1 = days[i-1], cumulative[i-1]
                    x2, y2 = days[i], cumulative[i]
                    if y2 != y1:
                        return x1 + (target - y1) * (x2 - x1) / (y2 - y1)
                    else:
                        return x1
        return days[-1]

def create_germination_curves(data, treatment_name="Treatment", save_path=""):
    plt.style.use('dark_background')
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
    if 'Day' in data.columns:
        days = data['Day']
    elif 'ti' in data.columns:
        days = data['ti']
    else:
        days = data.iloc[:, 0]
    rep_columns = [col for col in data.columns if col.startswith('R')]
    if not rep_columns:
        rep_columns = data.columns[1:]
    colors = ['#00b894', '#00cec9', '#0984e3', '#fdcb6e', '#d63031', '#e17055']

    # 1. Individual Cumulative Curves
    cumulative_data = []
    for i, col in enumerate(rep_columns[:6]):
        if col in data.columns:
            daily_data = data[col].fillna(0)
            cumulative = np.cumsum(daily_data)
            cumulative_data.append(cumulative.values)
            ax1.plot(days, cumulative, marker='o', linewidth=2.5,
                     label=f'Rep {i+1}', color=colors[i % len(colors)], markersize=4, alpha=0.8)
    ax1.set_xlabel('Days', fontsize=12, fontweight='bold', color='white')
    ax1.set_ylabel('Cumulative Germinated Seeds', fontsize=12, fontweight='bold', color='white')
    ax1.set_title(f'Individual Cumulative Germination - {treatment_name}', fontsize=14, fontweight='bold', color='white')
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # 2. Mean Cumulative Curve with Error Bars
    if cumulative_data:
        cumulative_array = np.array(cumulative_data).T
        mean_cumulative = np.mean(cumulative_array, axis=1)
        std_cumulative = np.std(cumulative_array, axis=1)
        ax2.plot(days, mean_cumulative, marker='o', linewidth=3, color='#00b894',
                 label='Mean', markersize=6)
        ax2.fill_between(days, mean_cumulative - std_cumulative,
                         mean_cumulative + std_cumulative, alpha=0.3, color='#00cec9',
                         label='±1 SD')
    ax2.set_xlabel('Days', fontsize=12, fontweight='bold', color='white')
    ax2.set_ylabel('Mean Cumulative', fontsize=12, fontweight='bold', color='white')
    ax2.set_title(f'Mean Cumulative Germination (±SD) - {treatment_name}', fontsize=14, fontweight='bold', color='white')
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    # 3. Individual Daily Germination
    for i, col in enumerate(rep_columns[:6]):
        if col in data.columns:
            daily_data = data[col].fillna(0)
            ax3.plot(days, daily_data, marker='s', linewidth=2,
                     label=f'Rep {i+1}', color=colors[i % len(colors)], markersize=4, alpha=0.8)
    ax3.set_xlabel('Days', fontsize=12, fontweight='bold', color='white')
    ax3.set_ylabel('Seeds Germinated per Day', fontsize=12, fontweight='bold', color='white')
    ax3.set_title(f'Daily Germination Frequency - {treatment_name}', fontsize=14, fontweight='bold', color='white')
    ax3.legend()
    ax3.grid(True, alpha=0.3)

    # 4. Mean Daily Germination with Error Bars
    rep_data = data[rep_columns].fillna(0)
    mean_daily = rep_data.mean(axis=1)
    std_daily = rep_data.std(axis=1)
    ax4.bar(days, mean_daily, alpha=0.7, color='#00b894', edgecolor='#00cec9',
            linewidth=1.5, label='Mean Daily Germination')
    ax4.errorbar(days, mean_daily, yerr=std_daily, fmt='none', color='#d63031',
                 capsize=4, capthick=2, label='±1 SD')
    ax4.set_xlabel('Days', fontsize=12, fontweight='bold', color='white')
    ax4.set_ylabel('Mean Seeds Germinated per Day', fontsize=12, fontweight='bold', color='white')
    ax4.set_title(f'Mean Daily Germination Pattern (±SD) - {treatment_name}', fontsize=14, fontweight='bold', color='white')
    ax4.legend()
    ax4.grid(True, alpha=0.3)

    for ax in [ax1, ax2, ax3, ax4]:
        ax.tick_params(axis='x', colors='white')
        ax.tick_params(axis='y', colors='white')
        for spine in ax.spines.values():
            spine.set_color('white')

    plt.tight_layout(pad=3.0)
    if save_path:
        filename = os.path.join(save_path, f'GermiTrack_Germination_Curves_{treatment_name}.png')
    else:
        filename = f'GermiTrack_Germination_Curves_{treatment_name}.png'
    plt.savefig(filename, dpi=300, bbox_inches='tight', facecolor='black')
    plt.savefig(filename.replace('.png', '.pdf'), bbox_inches='tight', facecolor='black')

    save_germination_curves_data(data, treatment_name, save_path)

    return fig, filename

# ==========================================
# 2. Modern Graphical Interface (Tkinter)
# ==========================================

class GermiTrackDesktop:
    def __init__(self):
        self.root = tk.Tk()
        self.root.title("GermiTrack v.1.4 - Pontes et al. (2025)")
        self.root.geometry("1050x750")
        self.root.configure(bg='#18191a')

        # Variables
        self.data = {}
        self.results = {}
        self.total_seeds = tk.IntVar(value=25)
        self.output_folder = tk.StringVar(value=str(Path.home() / "Desktop" / "GermiTrack_Results"))

        self.setup_gui()
        self.center_window()

    def center_window(self):
        self.root.update_idletasks()
        x = (self.root.winfo_screenwidth() // 2) - (self.root.winfo_width() // 2)
        y = (self.root.winfo_screenheight() // 2) - (self.root.winfo_height() // 2)
        self.root.geometry(f"+{x}+{y}")

    def setup_gui(self):
        style = ttk.Style()
        style.theme_use('clam')
        style.configure('TFrame', background='#18191a')
        style.configure('TLabel', background='#18191a', foreground='white', font=('Arial', 11))
        style.configure('TButton', font=('Arial', 11, 'bold'), padding=6)
        style.configure('TLabelframe', background='#23272e', foreground='white', font=('Arial', 12, 'bold'))
        style.configure('TLabelframe.Label', background='#23272e', foreground='white', font=('Arial', 12, 'bold'))

        main_frame = ttk.Frame(self.root, padding="24", style='TFrame')
        main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)
        main_frame.columnconfigure(1, weight=1)

        # Title and authors
        title_label = ttk.Label(main_frame, text="🌱 GermiTrack: Smart Analysis for Germination Data",
                               font=('Arial', 20, 'bold'), foreground='#00b894', background='#18191a')
        title_label.grid(row=0, column=1, columnspan=2, pady=(0, 10), sticky=tk.W)

        authors_label = ttk.Label(main_frame, text="Authors: Pontes et al. (2025)",
                                  font=('Arial', 11, 'italic'), foreground='#b2bec3', background='#18191a')
        authors_label.grid(row=1, column=1, columnspan=2, pady=(0, 18), sticky=tk.W)

        # Institutional note
        note_label = ttk.Label(main_frame, text="Developed by GOF/UFMS",
                               font=('Arial', 10, 'italic'), foreground='#636e72', background='#18191a')
        note_label.grid(row=2, column=0, columnspan=3, pady=(0, 15), sticky=tk.W)

        # File upload
        upload_frame = ttk.Labelframe(main_frame, text="Data Input", padding="12", style='TLabelframe')
        upload_frame.grid(row=3, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=(0, 12))
        upload_frame.columnconfigure(1, weight=1)
        ttk.Button(upload_frame, text="📁 Select File",
                   command=self.upload_file, width=22).grid(row=0, column=0, padx=(0, 10))
        self.file_label = ttk.Label(upload_frame, text="No file selected", foreground="#636e72", background='#23272e')
        self.file_label.grid(row=0, column=1, sticky=tk.W)

        # Parameters
        params_frame = ttk.Labelframe(main_frame, text="Parameters", padding="12", style='TLabelframe')
        params_frame.grid(row=4, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=(0, 12))
        ttk.Label(params_frame, text="Total Seeds per Replicate:", background='#23272e', foreground='white').grid(row=0, column=0, sticky=tk.W)
        ttk.Spinbox(params_frame, from_=1, to=1000, textvariable=self.total_seeds,
                    width=10).grid(row=0, column=1, sticky=tk.W, padx=(10, 0))
        ttk.Label(params_frame, text="Output Folder:", background='#23272e', foreground='white').grid(row=1, column=0, sticky=tk.W, pady=(10, 0))
        ttk.Entry(params_frame, textvariable=self.output_folder, width=52).grid(row=1, column=1,
                    sticky=(tk.W, tk.E), padx=(10, 5), pady=(10, 0))
        ttk.Button(params_frame, text="Browse",
                   command=self.select_output_folder).grid(row=1, column=2, pady=(10, 0))
        params_frame.columnconfigure(1, weight=1)

        # Action buttons
        buttons_frame = ttk.Frame(main_frame, style='TFrame')
        buttons_frame.grid(row=5, column=0, columnspan=3, pady=18)
        ttk.Button(buttons_frame, text="🔬 Analyze Data",
                   command=self.analyze_data, width=22).grid(row=0, column=0, padx=6)
        ttk.Button(buttons_frame, text="📊 Open Results Folder",
                   command=self.open_results_folder, width=22).grid(row=0, column=1, padx=6)

        # Progress bar
        self.progress = ttk.Progressbar(main_frame, mode='indeterminate')
        self.progress.grid(row=6, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=(10, 0))

        # Output log area
        output_frame = ttk.Labelframe(main_frame, text="Output Log", padding="12", style='TLabelframe')
        output_frame.grid(row=7, column=0, columnspan=3, sticky=(tk.W, tk.E, tk.N, tk.S), pady=(10, 0))
        output_frame.columnconfigure(0, weight=1)
        output_frame.rowconfigure(0, weight=1)
        self.output_text = scrolledtext.ScrolledText(output_frame, height=15, width=90, font=("Consolas", 10), bg='#23272e', fg='white', insertbackground='white')
        self.output_text.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        main_frame.rowconfigure(7, weight=1)

    def log_message(self, message):
        self.output_text.insert(tk.END, f"{message}\n")
        self.output_text.see(tk.END)
        self.root.update()

    def upload_file(self):
        filename = filedialog.askopenfilename(
            title="Select germination data file",
            filetypes=[
                ("Excel files", "*.xlsx *.xls"),
                ("CSV files", "*.csv"),
                ("Text files", "*.txt"),
                ("All files", "*.*")
            ]
        )
        if filename:
            try:
                self.log_message(f"Loading file: {os.path.basename(filename)}")
                if filename.endswith(('.xlsx', '.xls')):
                    excel_file = pd.ExcelFile(filename)
                    self.data = {}
                    for sheet in excel_file.sheet_names:
                        self.data[sheet] = pd.read_excel(filename, sheet_name=sheet)
                    self.log_message(f"✅ Excel file loaded!")
                    self.log_message(f"Sheets: {', '.join(excel_file.sheet_names)}")
                elif filename.endswith('.csv'):
                    self.data = {'Data': pd.read_csv(filename)}
                    self.log_message("✅ CSV file loaded!")
                elif filename.endswith('.txt'):
                    self.data = {'Data': pd.read_csv(filename, sep='\t')}
                    self.log_message("✅ TXT file loaded!")
                self.file_label.config(text=os.path.basename(filename), foreground="#00b894", background='#23272e')
                for sheet_name, df in self.data.items():
                    self.log_message(f"\n📊 Preview - {sheet_name}:")
                    self.log_message(f" Rows: {len(df)}, Columns: {len(df.columns)}")
                    self.log_message(f" Germination events: {df.iloc[:, 1:].sum().sum()}")
            except Exception as e:
                messagebox.showerror("Error", f"Error loading file:\n{str(e)}")
                self.log_message(f"❌ Error loading file: {str(e)}")

    def select_output_folder(self):
        folder = filedialog.askdirectory(title="Select output folder")
        if folder:
            self.output_folder.set(folder)

    def analyze_data(self):
        if not self.data:
            messagebox.showwarning("Warning", "Please upload a data file first!")
            return
        try:
            output_dir = Path(self.output_folder.get())
            output_dir.mkdir(parents=True, exist_ok=True)
            self.progress.start()
            self.log_message("\n🔬 Starting germination analysis...")
            all_results = {}
            for treatment, data in self.data.items():
                if len(data) > 0 and not data.empty:
                    self.log_message(f"\n📊 Analyzing treatment: {treatment}")
                    analyzer = GermiTrackAnalyzer(data, self.total_seeds.get())
                    results = analyzer.calculate_all_germination_parameters(data)
                    all_results[treatment] = results
                    self.log_message(f" Generating charts and exporting data...")
                    fig, plot_filename = create_germination_curves(data, treatment, str(output_dir))
                    plt.close(fig)
                    self.log_message(f" ✅ Charts and data saved: {os.path.basename(plot_filename)}")
            if all_results:
                self.log_message("\n💾 Saving complete Excel report...")
                excel_filename = output_dir / 'GermiTrack_Complete_Analysis.xlsx'
                with pd.ExcelWriter(excel_filename, engine='openpyxl') as writer:
                    for treatment, results in all_results.items():
                        results.to_excel(writer, sheet_name=f'{treatment}_Results', index=False)
                self.log_message("\n✅ Analysis completed successfully!")
                self.log_message(f"📁 Results saved in: {output_dir}")
                messagebox.showinfo("Analysis Complete",
                                    f"Analysis completed successfully!\n\nResults saved in:\n{output_dir}")
            self.progress.stop()
        except Exception as e:
            self.progress.stop()
            error_msg = f"Error during analysis: {str(e)}"
            self.log_message(f"❌ {error_msg}")
            messagebox.showerror("Error", error_msg)

    def open_results_folder(self):
        output_dir = Path(self.output_folder.get())
        if output_dir.exists():
            if sys.platform == "win32":
                os.startfile(output_dir)
            elif sys.platform == "darwin":
                os.system(f"open '{output_dir}'")
            else:
                os.system(f"xdg-open '{output_dir}'")
        else:
            messagebox.showwarning("Warning", "Results folder does not exist yet!")

    def run(self):
        self.root.mainloop()

def main():
    try:
        app = GermiTrackDesktop()
        app.run()
    except Exception as e:
        print(f"Error starting GermiTrack: {e}")
        input("Press Enter to exit...")

if __name__ == "__main__":
    main()
