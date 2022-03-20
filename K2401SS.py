# Keithly 2401 Solar Simulator class -- RLH 12/24/2021
"""
Run Solar cell Current-Voltage scans on Keithley Sourcemeter.
Also control Newport LSH-7320 lamp trigger input via a Labjack U3 digital output.

class K2401 functions:
    constructor: Connect to Keithley SMU on GPIB channel 24. 
    dll_dir = r'C:\WINDOWS\system32\ visa64.dll' by default.
    Other backends, such as PyVISA-py are possible.
    
    IVsweep(startv, stopv, stepv, rate='medium') : perform a voltage
    sweep while measuring current.  Rates of 'slow', 'medium' and 'fast'
    are suppored. They correspond to 1.0, 0.1, and 0.0 line voltage cycle
    delay between each data point.

    dark_light_IVsweep(startv, stopv, stepv, rate='medium', 
    area=1.0, plottitle='', verbose=1): Call IVsweep twice, once with
    the lamp off, and a second time with it on. The cell area is in units
    of cm^2. Use it to calculate J (mA/cm^2).

    save_data_and_plot(df, fig, path = 'C:\\Users\Public', 
    samplename = 'testymctest', cellnumber = '0', timenowstr = '', verbose='yes'):
    Save the data (pandas DataFrame) in a csv file and also save a plot of the dark
    and light IV sweeps.

    calc_and_save_parameters(df,  path = 'C:\\Users\Public', 
    samplename = 'testymctest', cellnumber = '0', timenowstr = '', 
    datetodaystr = '', saveparameters = 'yes', verbose = 1, timeseries = False,
    base_t = 0): Calculate the photovoltaic cell parameters -- Voc,
    Isc, Rsh (shunt resistance), Rser (series resistance), FF, PCE.
    Note: The PCE is calculated assuming the calibration of the LSH-7320,
    which is 100 mW/cm^2 between 400 and 1100 nm. AM1.5 100 mW/cm^2 has
    only 80.5 mW/cm^2 between 400-1100 nm. This factor is applied and no
    other correction factors are used. Cell parameters are saved to a csv file, 
    which is appended-to after subsequent measurement. If timeseries==True,
    a special timeseries_parameters file is created with an extra column,
    'Elapsed Time', which records the time of each measurement relative to
    the time that the '/data' folder was created.

 class LabJackU3:   
    Control the lamp through the LabJack FIO4 output set in 
    'digital output mode'. The high level is 3.5V, which 
    corresponds to the LSH-7320 being off. 
"""
class K2401:
    import matplotlib.pyplot as plt
    from scipy import stats

    def __init__(self, dll_dir=r'C:\WINDOWS\system32\visa64.dll'):
        """
        Initialize the device.
        """
        self.ctrl = 0 # The sourcemeter device.
        import os
        import sys
        import pyvisa as visa # PyVISA library
        print("K2401 setup:\n")
        print("Default Directory: ", os.getcwd())

        # Parts of this code block were derived from:
        # https://rfmw.em.keysight.com//DigitalPhotonics/flexdca/PG/Content/Topics/Python/a_python.htm
        rm = visa.ResourceManager(dll_dir)
        #rm = visa.ResourceManager('@py') # PyVISA-py backend
        print('\nVISA library version:\n  ', rm)
        print('\nPython version:\n  ', sys.version)
        print('\nList of instruments discovered:')
        # Loop over all resources and print the instrument names.
        i=1
        for key, value in rm.list_resources_info().items():
            print('\nInstrument ', str(i), ': ', key)
            print('  Interface type: ', value.interface_type)
            print('  Interface board number: ', value.interface_board_number)
            print('  Resource class: ', value.resource_class)
            print('  Resource name: ', value.resource_name)
            print('  Resource alias: ', value.alias)
            i += 1
        # Open the device.
        ADDRESS = 'GPIB1::25::INSTR'
        self.ctrl =  rm.open_resource(ADDRESS)
        self.ctrl.write('*IDN?')
        print("\n", self.ctrl.read())
        return 

    def IVsweep(self, startv, stopv, stepv, rate="medium"):
        from numpy import arange, ceil
        from pandas import DataFrame
        """
        Sweep the voltage, measure current.
        Do one dark sweep with LSH-7320 off, then one light.
        Return a dataframe with the Voltage and Current data for both dark and light
        """
        # self.ctrl IV Sweep for 2400 SourceMeter
        # See github.com/pmasi/self.ctrl-Python-Interface

        # Set delay time
        if (rate == "fast"):
            delay = "0.0"
        elif (rate == "slow"):
            delay = "1.0"
        else:
            delay = "0.1"

        # Variable assignment
        filename = "test_test_1"
        startvprime = float(startv)
        stopvprime = float(stopv)
        stepvprime = float(stepv)
        steps = ceil((stopv - startv) / stepv)
        startv_str = str(startv)
        stopv_str = str(stopv)
        stepv_str = str(stepv)

        self.ctrl.write("*RST")
        self.ctrl.timeout = 25000

        # Turn off concurrent functions and set sensor to current with fixed voltage
        self.ctrl.write(":SENS:FUNC:CONC OFF")
        self.ctrl.write(":SOUR:FUNC VOLT")
        self.ctrl.write(":SENS:FUNC 'CURR:DC' ")
        if (delay != "0.0"):
            self.ctrl.write(":SOUR:DEL ", delay)

        # Voltage starting, ending, and spacing values based on input
        self.ctrl.write(":SOUR:VOLT:STAR ", startv_str)
        self.ctrl.write(":SOUR:VOLT:STOP ", stopv_str)
        self.ctrl.write(":SOUR:VOLT:STEP ", stepv_str)
        self.ctrl.write(":SOUR:SWE:RANG AUTO")

        # Set compliance current (in A), sweep direction, and data acquisition
        self.ctrl.write(":SENS:CURR:PROT 1")
        self.ctrl.write(":SOUR:SWE:SPAC LIN")
        self.ctrl.write(":SOUR:SWE:POIN ", str(int(steps)))
        self.ctrl.write(":SOUR:SWE:DIR UP")
        self.ctrl.write(":TRIG:COUN ", str(int(steps)))
        self.ctrl.write(":FORM:ELEM CURR")

        # Set sweep mode and turn output on
        self.ctrl.write(":SOUR:VOLT:MODE SWE")
        self.ctrl.write(":OUTP ON")

        # Initiate sweep, collect ACSII current values, and turn output off
        result = self.ctrl.query(":READ?")
        yvalues = self.ctrl.query_ascii_values(":FETC?")
        self.ctrl.write(":OUTP OFF")
        self.ctrl.write(":SOUR:VOLT 0")

        # Create xvalues array.
        xvalues = arange(startv,stopv,stepv)

        return DataFrame({"V(V)": xvalues, "I(A)": yvalues})

    def dark_light_IVsweep(self, startv, stopv, stepv, rate="medium", 
    area=1.0, plottitle="", verbose=1):
        from numpy import arange
        from pandas import DataFrame
        import time
        import matplotlib.pyplot as plt
        """
        Sweep the voltage, measure current.
        Return a dataframe witht the Voltage and Current data 
        """
        if(verbose>0):
            print("LabJack U3 setup: ")
        U3 = LabjackU3(verbose=verbose)     # LabJack U3 device

        # Dark measurement
        if(verbose>0):
            print("Starting IV scans:")
        U3.set_lamp_state("off", verbose=verbose)
        start_t = time.time()
        dfd = self.IVsweep(startv, stopv, stepv, rate=rate)
        stop_t = time.time()
        if(verbose>0):
            print("Sweep rate: {:3.1f} mV/sec".format(1000*(stopv-startv)/(stop_t-start_t)))
        xvalues_dark = dfd["V(V)"]
        yvalues_dark = dfd["I(A)"]

        # Light measurement
        U3.set_lamp_state("on", verbose=verbose)
        time.sleep(1) # wait a few seconds for the lamp to stabilize
        start_t = time.time()
        dfl = self.IVsweep(startv, stopv, stepv, rate=rate)
        stop_t = time.time()
        if(verbose>0):
            print("Sweep rate: {:3.1f} mV/sec".format(1000*(stopv-startv)/(stop_t-start_t)))
        xvalues_light = dfl["V(V)"]
        yvalues_light = dfl["I(A)"]

        # Turn the lamp off
        U3.set_lamp_state("off", verbose=verbose)
        if(verbose>0):
            print("Data collection is done!")

        # Make the plot
        fig = plt.figure()
        plt.plot(xvalues_dark,1000*yvalues_dark/area, label="dark")
        plt.plot(xvalues_light,1000*yvalues_light/area, label="light")
        plt.xlabel(' Voltage (V)')
        plt.ylabel(' J (mA/cm$^2$)')
        plt.grid()
        plt.title(plottitle)
        plt.legend()
        if(verbose>0):
            plt.show()
        else:
            plt.close()

        # Make a new Dataframe with both Dark and Light data
        d = {plottitle + "_V_D": xvalues_dark, plottitle + "_J_D": yvalues_dark/area, 
             plottitle + "_V_L": xvalues_light, plottitle + "_J_L": yvalues_light/area }
        return fig, DataFrame(data=d)

    def save_data_and_plot(self, df, fig, path = "C:\\Users\Public", 
    samplename = "testymctest", cellnumber = "0", timenowstr = "", verbose="yes"):
        """Save data from the dataframe to a file. Show and save the plot from fig."""
        from datetime import datetime
        import matplotlib.pyplot as plt
        from scipy.interpolate import interp1d
        import time
        if (timenowstr == ""):
            timenow = datetime.now()
            timenowstr = timenow.strftime("%H%M%S")
        sample_cellnum_time_str = samplename + "_" + cellnumber  + "_" + timenowstr
        filename = path + "\\" + sample_cellnum_time_str
        df.to_csv(filename)
        print("Data saved to: " + filename + ".csv")
        fig.savefig(filename + ".png",format="png",
                bbox_inches ="tight",
                pad_inches = 1,
                transparent = True,
                facecolor ="w",
                edgecolor ='b',
                orientation ='landscape')
        plt.show()
        print("Plot saved to: " + filename + ".png")

    def calc_and_save_parameters(self, df,  path = "C:\\Users\Public", 
    samplename = "testymctest", cellnumber = "0", timenowstr = "", 
    datetodaystr = "", saveparameters = "yes", verbose = 1, timeseries = False,
    base_t = 0, I_ph = 1.0):
        """
        Calculate PV parameters and save to a file. Append if the file already exists.
        Note: AM1.5 has 80.5 mW/cm2 between 400-1100 nm. LSH-7320 calibration is 100 mW/cm2 at 1 sun.
        """
        import os
        import numpy as np
        import pandas as pd
        from scipy.interpolate import interp1d
        from datetime import datetime
        import time
        #* Calculate V_OC, J_SC, R_Sh, R_Ser, FF, PCE and save them to a common "parameters" file.
        v = df[samplename+"_"+cellnumber + "_V_L"]
        j = df[samplename+"_"+cellnumber + "_J_L"]
        #
        fj = interp1d(v,j,kind='cubic')
        j_sc = -1000.0*fj(0.0)
        vstep = 0.01
        r_sh = vstep/(fj(vstep)-fj(0.0))
        try:
            fv = interp1d(j,v,kind='cubic')
            v_oc = fv(0.0)
            jstep = 1e-6
            r_ser = (fv(jstep)-fv(0.0))/jstep
            pwr = 1000*(-v*j)
            pmax = np.max(pwr)
            pmaxpos = np.where(pwr == pmax)
            mpp = v[pmaxpos[0][0]] 
            ff = pmax/(v_oc*j_sc)
            pce = 0.805*pmax/(100.0*I_ph) # AM1.5 has 80.5 mW/cm2 between 400-1100 nm. LSH-7320 calibration is 100 mW/cm2 at 1 sun.
        except:
            print("Warning: interp1d error while calculating v_oc")
            v_oc = 0.0 # not the correct value. This is just to handle the error.
            r_ser = 0.0
            pmax = 0.0
            mpp = 0.0
            ff = 0.0
            pce = 0.0
        # 
        if (verbose == 1):
            print("I_ph = {:5.3f} V".format(I_ph))
            print("V_oc = {:5.3f} V".format(v_oc))
            print("mpp = {:5.3f} V".format(mpp))
            print("J_sc = {:5.3f} mA".format(j_sc))
            print("R_sh = {:5.1f} Ω-cm^2".format(r_sh))
            print("R_ser = {:5.2f} Ω-cm^2".format(r_ser))
            print("FF = {:5.3f}".format(ff))
            print("PCE = {:5.2f}%".format(100*pce))
        #
        #### Save parameters
        sample_cellnum_time_str = samplename +"_"+ cellnumber + "_"+ timenowstr
        if (timenowstr == ""):
            timenow = datetime.now()
            timenowstr = timenow.strftime("%H%M%S")
        if (timeseries == True):
            filename = path + "\\" + "timeseries_parameters_" + datetodaystr
            data = {"Elapsed Time": (time.time() - base_t)/60.0,
                    "Sample_Cell#":  sample_cellnum_time_str,
                    "I_ph(suns)": [I_ph],
                    "V_oc(V)":  [v_oc],
                    "mpp(V)":  [mpp],
                    "J_sc(mA)":  [j_sc],
                    "R_sh":  [r_sh],
                    "R_ser": [r_ser],
                    "FF":    [ff],
                    "PCE(%)":   [100*pce]}
        else:
            filename = path + "\\" + "parameters_" + datetodaystr
            data = {"Sample_Cell#":  sample_cellnum_time_str,
                    "I_ph(suns)": [I_ph],
                    "V_oc(V)":  [v_oc],
                    "mpp(V)":  [mpp],
                    "J_sc(mA)":  [j_sc],
                    "R_sh":  [r_sh],
                    "R_ser": [r_ser],
                    "FF":    [ff],
                    "PCE(%)":   [100*pce]}
        #
        par_df_new = pd.DataFrame(data)
        if (saveparameters != 'no'):
            # Check whether the specified file exists or not
            isExist = os.path.exists(filename)
            if isExist:
                # Open the file and read the dataframe to recover previously saved parameters.
                # Then add the new parameters and save the file again, overwriting the old version.
                par_df = pd.read_csv(filename, index_col=False)
                par_df = pd.concat([par_df, par_df_new], ignore_index = True, axis = 0)
                par_df.to_csv(filename, index=False)
                print("Appended parameters to:", filename)
            else:
                # Create a new dataframe with the right headings.
                par_df = par_df_new
                par_df.to_csv(filename, index=False)
                print("Created: ", filename)
        else:
            print("Parameters NOT saved.")
        #
        return par_df_new



class LabjackU3:

    def __init__(self, verbose=1):
        import u3
        self.d = u3.U3()
        val = self.d.configU3(FIOAnalog = 0)
        if(verbose>0):
            print(val) # digital mode
        return

    def set_lamp_state(self, mystate, verbose=1):
        # turn LSH-7320 lamp on or off
        off = 1; on = 0
        #
        if (mystate == "on") :
            newstate = on 
            textns = "on"
        else: 
            newstate = off
            textns = "off"
        #
        self.d.setFIOState(4, state = newstate)
        #
        if self.d.getDIOState(4) == newstate:
            if(verbose>0):
                print("Lamp", textns)
        else: 
            print("Warning: lamp didn't switch properly")