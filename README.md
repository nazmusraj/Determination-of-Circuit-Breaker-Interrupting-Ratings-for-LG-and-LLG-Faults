# Power System Load Flow and Fault Analysis Tool

A comprehensive MATLAB-based software tool for steady-state load flow and unsymmetrical fault analysis of power systems. This project was developed as part of the EEE 306: Power System I Laboratory course at the Bangladesh University of Engineering and Technology (BUET).

## Project Overview

This repository contains a collection of MATLAB scripts and a Graphical User Interface (GUI) designed to perform two critical power system studies:

1.  **Load Flow Analysis:** Determines the steady-state operating conditions of a power system, including bus voltages, phase angles, and power flows. The tool includes implementations of three classical methods:
    * Fast-Decoupled Load Flow (FDLF)
    * Newton-Raphson (N-R)
    * Gauss-Seidel (G-S)

2.  **Unsymmetrical Fault Analysis:** Using the pre-fault conditions from the load flow, the tool calculates fault currents for the two most common types of unbalanced faults:
    * Single Line-to-Ground (LG) Fault
    * Double Line-to-Ground (LLG) Fault

The primary outcome is the determination of the **Symmetrical Interrupting MVA Rating** for circuit breakers, a crucial parameter for designing robust protection schemes.

## Features

- **Multiple Load Flow Solvers:** Compare results from FDLF, Newton-Raphson, and Gauss-Seidel methods.
- **Interactive Fault Analysis:** After a successful load flow, users can interactively apply LG or LLG faults to any bus in the system.
- **GUI for Ease of Use:** A user-friendly Graphical User Interface (GUI) built with MATLAB App Designer allows for easy data loading and analysis without touching the code.
- **Excel Data Input:** System data (bus and line parameters) is managed in simple `.xlsx` files, making it easy to test different systems.
- **Detailed Results:** The tool provides comprehensive output for both load flow (bus data, losses) and fault analysis (sequence currents, phase currents, breaker ratings).
- **Export Functionality:** Load flow results can be exported to a new Excel file directly from the scripts or the GUI.

## Repository Contents

-   `FDLF_Analysis.m`: MATLAB script for FDLF Load Flow and Fault Analysis.
-   `Newton.m`: MATLAB script for Newton-Raphson Load Flow and Fault Analysis.
-   `Gauss.m`: MATLAB script for Gauss-Seidel Load Flow and Fault Analysis.
-   `PowerSystemAnalyzerGUI.mlapp`: The MATLAB App Designer file for the graphical user interface.
-   `/Project_Report.pdf` final project report
-   `/GUI1/`: Matlab GUI first version. Inside all necessary files are given.
-   `/GUI2/`: Matlab GUI Second version. Inside all necessary files are given.
-   `/GUI3/`: Matlab GUI third version. Inside all necessary files are given.
-   `README.md`: This file.

## How to Use

### Prerequisites

-   MATLAB (R2021a or newer recommended)
-   MATLAB App Designer (for the GUI)

### 1. Using the MATLAB Scripts (e.g., `FDLF_Analysis.m`)

1.  Ensure the script and the `/Data_Files/` folder are in the same directory.
2.  Open the desired script (e.g., `FDLF_Analysis.m`) in MATLAB.
3.  In the "Data Input" section, select the system you want to analyze by uncommenting the appropriate `fileName`.
    ```matlab
    % Select the system to analyze
    fileName = 'IEEE_5.xlsx';
    % fileName = 'IEEE_9.xlsx';
    ```
4.  Run the script. The load flow results will be displayed in the command window and saved to a new Excel file.
5.  The script will then enter the interactive fault analysis mode. Follow the prompts in the command window to test fault scenarios.

### 2. Using the GUI (`PowerSystemAnalyzerGUI.mlapp`)

1.  Open MATLAB.
2.  Navigate to the project directory.
3.  Double-click `PowerSystemAnalyzerGUI.mlapp` to open it in App Designer, then click the "Run" button.
4.  **Step 1:** Click the "Browse..." button to select an Excel data file (e.g., `IEEE_9.xlsx`).
5.  **Step 2:** Click the "Run Load Flow & Prepare" button. The load flow results will appear in the table.
6.  **Step 3:** Navigate to the "Fault Analysis" tab. Select the desired Fault Bus, Fault Type, and Fault Impedance.
7.  **Step 4:** Click the "Calculate Fault" button to see the results.

## Data Format

The input Excel files must contain two sheets named `BusData` and `LineData` with the following column structure:

**`BusData` Sheet:**
| Bus No. | Type | Vmag (p.u.) | Angle (deg) | PL (MW) | QL (MVAr) | PG (MW) | QG (MVAr) |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |

-   **Type:** 1 for Slack, 2 for PV, 3 for PQ.

**`LineData` Sheet:**
| From Bus | To Bus | R (p.u.) | X (p.u.) | B/2 (p.u.) | Tap |
| :--- | :--- | :--- | :--- | :--- | :--- |

## Contributors

This project was developed by:

-   **Md. Nazmus Sakib** (2106061) - *Team Lead, Main Coding & GUI*
-   & Others...

## License

This project is licensed under the MIT License. See the `LICENSE` file for details.
