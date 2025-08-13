GermiTrack Desktop v1.4
Smart Analysis for Germination Data

Overview
GermiTrack is a professional desktop application for the statistical analysis of seed germination experiments. It automates the calculation of all major germination parameters, generates publication-ready charts, and exports comprehensive results in Excel and CSV formats. Designed for researchers, agronomists, and seed technologists, GermiTrack delivers robust, reproducible, and transparent germination data analysis.

Authors & Institutional Note
Authors: Montcharles da Silva Pontes, Michele Aparecida dos Santos Nobrega, Samuel Leite de Oliveira, Anderson Rodrigues Lima Caires, 2025

Citation: Pontes et al. GermiTrack 1.4 - Smart Analysis for Germination Data, Campo Grande, MS, Brazil. GOF/UFMS, 2025.

Developed by: GOF/UFMS (Optics and Photonics Group, Federal University of Mato Grosso do Sul)

Features
Calculates all major germination parameters (germinability, mean germination time, synchrony, uncertainty, etc.)

Generates four main germination charts (individual/mean cumulative and daily curves, with standard deviation)

Exports all chart data and parameters to Excel and CSV

Modern, user-friendly graphical interface (Tkinter, dark mode)

Ready for publication and reporting


System Requirements
Operating System: Windows 10/11 (64-bit)

No Python required: All dependencies are bundled in the executable

Disk Space: At least 500 MB free for analyses and outputs

How to Use
1. Download and Extract
Copy the executable file GermiTrack_Desktop_1.4.exe to a folder of your choice.

(Optional) Place your germination data files (Excel, CSV, TXT) in the same or another convenient folder.

2. Run the Application
Double-click GermiTrack_Desktop_1.4.exe to start.

No installation or administrator rights are required.

3. Data Input
Click "Select File" and choose your data file.

Supported formats: .xlsx, .xls, .csv, .txt

Data should be organized with the first column as Day/Time and subsequent columns as replicates (R1, R2, ...).

4. Set Parameters
Enter the total number of seeds per replicate.

Choose or confirm the output folder for results.

5. Analyze
Click "Analyze Data".

The application will:

Calculate all germination parameters for each replicate and treatment.

Generate and save all charts (PNG and PDF).

Export all underlying data and results to Excel and CSV files.

Save a complete Excel report with all results.

6. Access Results
Click "Open Results Folder" to view all output files.

Output Files
Excel & CSV:

Individual and mean cumulative germination curves

Daily frequencies (individual and mean)

Standard deviations for all curves

Charts:

PNG and PDF files for each treatment

Comprehensive Excel Report:

All calculated parameters

Correlation matrix

Germination Parameters Calculated
Parameter	Description
Germinability (%)	Percentage of seeds that germinated out of the total sown.
Mean Germination Time (MGT)	Weighted average of germination times; lower values indicate faster germination.
Coefficient of Variation (CVt)	Relative dispersion of germination times; measures uniformity.
Mean Germination Rate (MGR)	Inverse of MGT; higher values indicate faster germination.
Uncertainty Index (U)	Shannon entropy-based index; higher values indicate greater temporal dispersion.
Synchrony Index (Z)	Degree of simultaneous germination; values near 1 mean highly synchronized germination.
Variance of Germination Time	Quantifies variability of germination times.
Standard Deviation	Dispersion of germination times around the mean.
Maguire Speed Index	Sum of (germinated seeds/day); higher values indicate rapid, concentrated germination.
T50 (Time to 50%)	Time required to reach 50% of final germination.
ArcSin Transformation	Statistical transformation of germinability for parametric analyses.
Troubleshooting
No file selected:
Make sure to upload a data file before running the analysis.

Permission Denied / Cannot Overwrite File:
Close any open output files before re-running analyses.

Antivirus Warning:
Some antivirus programs may flag new executables. Allow or whitelist the file if prompted.

Error: No module named 'matplotlib.backends.backend_pdf':
Use the provided .exe file; this error should not occur in the compiled version.

Contact
For questions, suggestions, or bug reports, contact:

Dr. Montcharles S. Pontes
E-mail subject: GermiTrack
E-mail: montcharles.pontes@gmail.com
GOF/UFMS (Optics and Photonics Group, Federal University of Mato Grosso do Sul)


License
This software is intended for academic and research use.
All rights reserved to the authors and GOF/UFMS.

Enjoy your professional germination data analysis with GermiTrack!
