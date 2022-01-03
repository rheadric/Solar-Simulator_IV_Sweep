# Solar-Simulator_IV

Control Keithley 2401 + LabJack U3 + Newport LSH-7320 solar simulator via python

K2401SS.py contains the code to communicate with the Keithley unit over GPIB, and with the Labjack over USB. The Labjack FIO4 is connected to the
trigger of the lamp to turn it on or off via a 3.5V TTL signal.

Solar_Simulator_IV_sweep_with_lamp_control_V2.ipynb: The IV-sweep Jupyter notebook implement Current-Voltage sweeps with the Keithley, while controlling whether the lamp is on or off for dark/light sweeps. A series of multiple sweeps can be produced, which is intended for PV cell degredation testing.

IV_sweep_Timescan_Monitor_V1.ipynb: A Monitor Jupyter notebook is also included with an animated plot that tracks PV PCE as a function of time, using the data generated by the first notebook.

The instructions file contains setup details for Windows and were tested on Windows 11.
