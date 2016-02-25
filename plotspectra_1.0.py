from matplotlib.figure import Figure 
import matplotlib
matplotlib.use('AGG')
matplotlib.rcParams['backend.qt4']='PySide'
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas 
from PySide import QtGui
import sys 
from numpy import *
import numpy as np
import matplotlib.cm as cm
import glob

class TabDialog(QtGui.QDialog): 
    def __init__(self): 
        super(TabDialog, self).__init__() 
        self.setGeometry(250, 200, 800, 700)
        self.setWindowTitle('PlotSpectra')
        self.setWindowIcon(QtGui.QIcon('logo.png'))  
        #############################  create the gas data loading layer
        gasLayout = QtGui.QHBoxLayout()
        gasLabel = QtGui.QLabel('Select Gas')
        self.gasCombo = QtGui.QComboBox(self)
        moleculeSpecies = molecules().moleculeSpecies()
        for mol in moleculeSpecies:
            self.gasCombo.addItem(mol)
        self.gasCombo.setMaxVisibleItems(50)
        loadButton = QtGui.QPushButton("Load") 
        loadButton.pressed.connect(self.LoadStatus)
        loadButton.clicked[bool].connect(self.LoadData) 
        gasLayout.addWidget(gasLabel)
        gasLayout.addWidget(self.gasCombo)
        gasLayout.addWidget(loadButton)
        self.loadBox = QtGui.QLabel('No gas selected')
        gasLayout.addWidget(self.loadBox)
        gasLayout.addStretch(1)       
        conditionLayout1 = QtGui.QHBoxLayout()
        conditionLayout2 = QtGui.QHBoxLayout()
        waveminLabel = QtGui.QLabel('Minimum Wavelength (um)')
        wavemaxLabel = QtGui.QLabel('Maximum Wavelength (um)')
        tempLabel = QtGui.QLabel('Temperature (oC)')
        pressureLabel = QtGui.QLabel('pressure (atm)')
        nptsLabel = QtGui.QLabel('number of points')
        concentrationLabel = QtGui.QLabel('Gas Percentage in Air')
        self.waveminBox = QtGui.QLineEdit('1')
        self.waveminBox.setFixedWidth(100)
        self.wavemaxBox = QtGui.QLineEdit('10')
        self.wavemaxBox.setFixedWidth(100)
        self.tempBox = QtGui.QLineEdit('27')
        self.tempBox.setFixedWidth(100)
        self.pressureBox = QtGui.QLineEdit('1')
        self.pressureBox.setFixedWidth(100)
        self.nptsBox = QtGui.QLineEdit('20000')
        self.nptsBox.setFixedWidth(100)
        self.concentrationBox = QtGui.QLineEdit('100')
        self.nptsBox.setFixedWidth(100)
        conditionLayout1.addWidget(waveminLabel)
        conditionLayout1.addWidget(self.waveminBox)
        conditionLayout1.addWidget(wavemaxLabel)
        conditionLayout1.addWidget(self.wavemaxBox)
        conditionLayout1.addWidget(nptsLabel)
        conditionLayout1.addWidget(self.nptsBox)
        conditionLayout1.addStretch(1)       
        conditionLayout2.addWidget(tempLabel)
        conditionLayout2.addWidget(self.tempBox)
        conditionLayout2.addWidget(pressureLabel)
        conditionLayout2.addWidget(self.pressureBox)
        conditionLayout2.addWidget(concentrationLabel)
        conditionLayout2.addWidget(self.concentrationBox)
        conditionLayout2.addStretch(1)
        #############################  create the range selection layer
        xlowLabel = QtGui.QLabel('XMin')
        xhighLabel = QtGui.QLabel('XMax')
        self.xlowBox = QtGui.QLineEdit()
        self.xlowBox.setFixedWidth(100)
        self.xhighBox = QtGui.QLineEdit()
        self.xhighBox.setFixedWidth(100)
        ylowLabel = QtGui.QLabel('yMin')
        yhighLabel = QtGui.QLabel('yMax')
        self.ylowBox = QtGui.QLineEdit()
        self.ylowBox.setFixedWidth(100)
        self.yhighBox = QtGui.QLineEdit()
        self.yhighBox.setFixedWidth(100)
        self.yautoscale = QtGui.QCheckBox('Autoscale')
        self.yautoscale.setChecked(True)
        rangeLayout = QtGui.QHBoxLayout()
        rangeLayout.addWidget(xlowLabel)
        rangeLayout.addWidget(self.xlowBox)
        rangeLayout.addWidget(xhighLabel)
        rangeLayout.addWidget(self.xhighBox)
        rangeLayout.addStretch(1)
        rangeLayout.addWidget(self.yautoscale)
        rangeLayout.addWidget(self.ylowBox)
        rangeLayout.addWidget(yhighLabel)
        rangeLayout.addWidget(self.yhighBox)
        ##############################    #create button layer
        buttonLayout = QtGui.QHBoxLayout() 
        self.plotButton = QtGui.QPushButton("Figure Plot") 
        self.plotButton.pressed.connect(self.PlotStatus) 
        self.plotButton.clicked.connect(self.Update) 
        buttonLayout.addWidget(self.plotButton) 
        self.plotButton.setDefault(True)        
        #############################  create the unit selection layer
        selectionLayout = QtGui.QHBoxLayout()
        unitSelection = QtGui.QButtonGroup(self)
        self.check_cm = QtGui.QRadioButton("cm-1", self)
        self.check_um = QtGui.QRadioButton("um", self)
        unitSelection.addButton(self.check_cm)
        unitSelection.addButton(self.check_um)
        selectionLayout.addWidget(self.check_cm)
        selectionLayout.addWidget(self.check_um)
        self.check_um.setChecked(True)

        ##############################    #create plot options layer
        seperateLayer1 = QtGui.QHBoxLayout()
        seperateLayer1.addWidget(QtGui.QLabel('----------------------------------------------'))
        yaxisLayout = QtGui.QHBoxLayout()
        yaxisSelection = QtGui.QButtonGroup(self)
        self.yaxis = QtGui.QComboBox(self)
        self.yaxis.addItem('Absorption')
        self.yaxis.addItem('Transmission')
        self.yaxis.addItem('Cross-section (cm2)')
        yaxisLayout.addWidget(self.yaxis)
        self.pathlengthLabel = QtGui.QLabel('Pathlength (cm)')
        self.pathlengthBox = QtGui.QLineEdit('10')
        self.pathlengthBox.setFixedWidth(100)
        yaxisLayout2 = QtGui.QHBoxLayout()
        yaxisLayout2.addWidget(self.pathlengthLabel)
        yaxisLayout2.addWidget(self.pathlengthBox)
        self.pathlengthLabel.setVisible(True)
        self.pathlengthBox.setVisible(True)
        self.yaxis.activated[str].connect(self.EnablePathLength)
        yaxisLayout2.addStretch(1) 
        #########################  put layout together
        totalLayout = QtGui.QHBoxLayout()
        sideLayout =  QtGui.QVBoxLayout()  
        sideLayout2 = QtGui.QVBoxLayout()
        separatorLabel1 = QtGui.QLabel('||')
#        separatorLabel2 = QtGui.QLabel('||')
#        separatorLabel3 = QtGui.QLabel('||')
#        separatorLabel4 = QtGui.QLabel('||')
#        separatorLabel5 = QtGui.QLabel('||')
        sideLayout2.addWidget(separatorLabel1)
#        sideLayout2.addWidget(separatorLabel2)
#        sideLayout2.addWidget(separatorLabel3)
#        sideLayout2.addWidget(separatorLabel4)
#        sideLayout2.addWidget(separatorLabel5)
        mainLayout = QtGui.QVBoxLayout()  
        mainLayout.addLayout(gasLayout)
        mainLayout.addWidget(AnalyzeTab()) 
        mainLayout.addLayout(conditionLayout1)
        mainLayout.addLayout(conditionLayout2)
        mainLayout.addLayout(rangeLayout) 
        mainLayout.addLayout(buttonLayout)
        sideLayout.addLayout(selectionLayout) 
        sideLayout.addLayout(seperateLayer1) 
        sideLayout.addLayout(yaxisLayout)
        sideLayout.addLayout(yaxisLayout2)
        sideLayout.addStretch(1) 
        totalLayout.addLayout(sideLayout)
        totalLayout.addLayout(sideLayout2)
        totalLayout.addLayout(mainLayout)
        self.setLayout(totalLayout) 
        self.currUnit = 1 # um = 1, cm = 0
        self.parameters = [self.check_um,self.xlowBox,self.xhighBox,self.currUnit,self.ylowBox,self.yhighBox, self.yaxis, self.gasCombo, self.pathlengthBox,self.tempBox,self.pressureBox, self.yautoscale] 
#==============================================================================
#     [0:self.check_um,1:self.xlowBox,2:self.xhighBox,3:currUnit,4:self.ylowBox,5:self.yhighBox, 
#        6:self.yaxis, 7:self.gasCombo, 8:self.pathlengthBox,9:self.tempBox,10:self.pressureBox,
#        11: self.yautoscale]
#       
#==============================================================================
    
    def EnablePathLength(self,input1):
        if input1 == 'Cross-section (cm2)':
            self.pathlengthLabel.setVisible(False)
            self.pathlengthBox.setVisible(False)  
        else:
            self.pathlengthLabel.setVisible(True)
            self.pathlengthBox.setVisible(True)    
    def PlotStatus(self):
        self.plotButton.setText('Plotting can take some time... Please wait...')
        self.plotButton.setStyleSheet("color: rgb(255, 0, 255);")
    def Update(self): 
        self.data = getData().calculate_hitran_xsec(self.data_raw, wavemin=float(self.waveminBox.text()), wavemax=float(self.wavemaxBox.text()), npts=int(self.nptsBox.text()), units='cm^2', temp=273+float(self.tempBox.text()), pressure=float(self.pressureBox.text()), concentration = self.concentrationBox)    
        canvas.update(self.data,self.parameters) 
        self.plotButton.setText("Figure Plot")
        self.plotButton.setStyleSheet("color: rgb(0, 0, 0);")
    def LoadStatus(self):
        self.loadBox.setText("loading data and plotting, please wait... \nsome data are huge so it could take tens of seconds")
        self.loadBox.setStyleSheet("color: rgb(255, 0, 255);")
    def LoadData(self):
        self.check_um.setChecked(True)
        self.parameters[3] = 1
        self.selectMol=self.gasCombo.currentText()
        self.data_raw = getData().read_hitran2012_parfile(molecules().filename(self.selectMol), self.loadBox)
        self.data = getData().calculate_hitran_xsec(self.data_raw, wavemin=float(self.waveminBox.text()), wavemax=float(self.wavemaxBox.text()), npts=int(self.nptsBox.text()), units='cm^2', temp=273+float(self.tempBox.text()), pressure=float(self.pressureBox.text()), concentration = self.concentrationBox)    
        canvas.FirstPlot(self.data, self.parameters)
        self.plotButton.setDefault(True)
        self.xlowBox.setText(str(min(canvas.ax.get_xlim())))
        self.xhighBox.setText(str(max(canvas.ax.get_xlim())))
        self.ylowBox.setText(str(min(canvas.ax.get_ylim())))
        self.yhighBox.setText(str(max(canvas.ax.get_ylim())))
        self.loadBox.setStyleSheet("color: rgb(0, 0, 0);")
        
class DataHolder(): 
    def __init__(self): 
        self.fig = Figure()
        self.ax = self.fig.add_subplot(111) 
        self.ax.set_title("Infrared Spectrum")
        self.ax.set_ylabel('Infrared Spectrum')
        self.ax.set_xlabel('Wavelength (um)')
        self.canvas = FigureCanvas(self.fig)
    def ProcessYaxis(self,ydata,parameters):
        density = float(parameters[10].text())*101325/1.38E-23/(float(parameters[9].text())+273)*1E-6 #measure in cm-3
        Trans = np.exp(-ydata*density*float(parameters[8].text()))
        yaxis = parameters[6].currentText()
        if yaxis == 'Cross-section (cm2)':
            return ydata
        elif yaxis == 'Transmission':
            return Trans
        else:
            return (1-Trans)
    def FirstPlot(self,data,parameters):
        self.ax.cla()
        self.yaxis = parameters[6].currentText()
        self.mol = parameters[7].currentText()
        self.ax.plot(data[0], self.ProcessYaxis(data[1],parameters), label=self.mol) 
        self.ax.legend(bbox_to_anchor=(1, 1))
        self.ax.set_title("Infrared Spectrum - "+self.yaxis)
        self.ax.set_ylabel(self.yaxis)
        self.ax.set_xlabel('Wavelength (um)')
        self.canvas.draw()
        self.ax.set_xlim([min(data[0]),max(data[0])])
    def update(self,data,parameters): 
        self.ax.cla()
        self.check_um = parameters[0]
        self.xmin = float(parameters[1].text())
        self.xmax = float(parameters[2].text())
        self.ymin = float(parameters[4].text())
        self.ymax = float(parameters[5].text())
        if self.check_um.isChecked():
            unit_um = True
            if parameters[3] == 0:
                self.xmin1 = self.xmin
                self.xmin = 10000/self.xmax
                self.xmax = 10000/self.xmin1
            parameters[3] = 1
            self.ax.set_xlabel('Wavelength (um)')
            self.ax.plot(data[0], self.ProcessYaxis(data[1],parameters), label=self.mol) 
        else:
            unit_um = False
            if parameters[3] == 1:
                self.xmin1 = self.xmin
                self.xmin = 10000/self.xmax
                self.xmax = 10000/self.xmin1
            parameters[3] = 0
            self.ax.set_xlabel('Wavenumber (cm-1)')        
            self.ax.plot(10000/data[0], self.ProcessYaxis(data[1],parameters), label=self.mol)       
        self.ax.legend(bbox_to_anchor=(1, 1))
        self.ax.set_ylabel(self.yaxis)
        self.ax.set_title("Infrared Spectrum - "+self.yaxis)        
        parameters[1].setText(str("%.3f"%(self.xmin)))
        parameters[2].setText(str("%.3f"%(self.xmax)))
        self.ax.set_xlim([self.xmin,self.xmax])
        if parameters[11].isChecked() == True:
            self.ax.autoscale(True)
        else:
            self.ax.set_ylim([self.ymin,self.ymax])
        self.canvas.draw()
        
class AnalyzeTab(QtGui.QWidget): 
    def __init__(self): 
        QtGui.QWidget.__init__(self) 
        analyzeLayout = QtGui.QGridLayout() 
        analyzeLayout.addWidget(canvas.canvas) 
        self.setLayout(analyzeLayout) 

class molecules():
    def __init__(self):
        self.moleculeDict = {'01':'H2O','02':'CO2','03':'O3','04':'N2O','05':'CO','06':'CH4','07':'O2','08':'NO',
              '09':'SO2','10':'NO2','11':'NH3','12':'HNO3','13':'OH','14':'HF','15':'HCl','16':'HBr',
             '17':'HI','18':'ClO','19':'OCS','20':'H2CO','21':'HOCl','22':'N2','23':'HCN','24':'CH3Cl',
             '25':'H2O2','26':'C2H2','27':'C2H6','28':'PH3','29':'COF2','30':'SF6','31':'H2S','32':'HCOOH',
             '33':'HO2','34':'O','35':'ClONO2','36':'NO+','37':'HOBr','38':'C2H4','39':'CH3OH','40':'CH3Br',
             '41':'CH3CN','42':'CF4','43':'C4H2','44':'HC3N','45':'H2','46':'CS','47':'SO3'}
    def filename(self,molecule):
        self.inv_moleculeDict = {v: k for k, v in self.moleculeDict.items()}
        return glob.glob('./par/' + self.inv_moleculeDict[molecule] +'*')[0]
    def moleculeSpecies(self):
        return sort(list(self.moleculeDict.values()))

        
        
class getData():
    def read_hitran2012_parfile(self, filename, loadbox):
        filehandle = open(filename, 'r')
        data = {'M':[],               ## molecule identification number
            'I':[],               ## isotope number
            'linecenter':[],      ## line center wavenumber (in cm^{-1})
            'S':[],               ## line strength, in cm^{-1} / (molecule m^{-2})
            'Acoeff':[],          ## Einstein A coefficient (in s^{-1})
            'gamma-air':[],       ## line HWHM for air-broadening
            'gamma-self':[],      ## line HWHM for self-emission-broadening
            'Epp':[],             ## energy of lower transition level (in cm^{-1})
            'N':[],               ## temperature-dependent exponent for "gamma-air"
            'delta':[],           ## air-pressure shift, in cm^{-1} / atm
            'Vp':[],              ## upper-state "global" quanta index
            'Vpp':[],             ## lower-state "global" quanta index
            'Qp':[],              ## upper-state "local" quanta index
            'Qpp':[],             ## lower-state "local" quanta index
            'Ierr':[],            ## uncertainty indices
            'Iref':[],            ## reference indices
            'flag':[],            ## flag
            'gp':[],              ## statistical weight of the upper state
            'gpp':[]}             ## statistical weight of the lower state
        for line in filehandle:
            if (len(line) < 160):
                loadbox.setText('The imported file ("' + filename + '") does not appear to be a HITRAN2012-format data file.')
            data['M'].append(uint(line[0:2]))
            data['I'].append(uint(line[2]))
            data['linecenter'].append(float64(line[3:15]))
            data['S'].append(float64(line[15:25]))
            data['Acoeff'].append(float64(line[25:35]))
            data['gamma-air'].append(float64(line[35:40]))
            data['gamma-self'].append(float64(line[40:45]))
            data['Epp'].append(float64(line[45:55]))
            data['N'].append(float64(line[55:59]))
            data['delta'].append(float64(line[59:67]))
            data['Vp'].append(line[67:82])
            data['Vpp'].append(line[82:97])
            data['Qp'].append(line[97:112])
            data['Qpp'].append(line[112:127])
            data['Ierr'].append(line[127:133])
            data['Iref'].append(line[133:145])
            data['flag'].append(line[145])
            data['gp'].append(line[146:153])
            data['gpp'].append(line[153:160])

        filehandle.close()
        loadbox.setText('Finished Loading' + filename + '\ndata loaded and plotted successfully')

        for key in data:
            data[key] = array(data[key])
        return(data)
        
    def calculate_hitran_xsec(self,data, wavemin=None, wavemax=None, npts=20001, units='cm^2', temp=296.0, pressure=1.0, concentration=100):
#==============================================================================
#     '''
#     Given the HITRAN data (line centers and line strengths) for a molecule, digitize the result into a spectrum of
#     absorption cross-section in units of cm^2.
#     Parameters
#     ----------
#     data : dict of ndarrays
#         The HITRAN data corresponding to a given molecule.
#     wavemin : float, optional
#         The minimum wavelength os the spectral region of interest.
#     wavemax : float, optional
#         The maximum wavelength os the spectral region of interest.
#     units : str, optional
#         A string describing in what units of the output cross-section should be given in. Choices available are:
#         {'cm^2/mole', 'cm^2.ppm', 'm^2/mole', 'm^2.ppm', 'm^2', cm^2}.
#     temp : float
#         The temperature of the gas, in Kelvin.
#     pressure : float
#         The pressure of the gas, in atmospheres.
#     Returns
#     -------
#     waves : ndarray
#         The wavelengths at which the cross-section data is evaluated.
#     xsec : array_like
#         The mean absorption cross-section (in cm^2) per molecule, evaluated at the wavelengths given by input `waves`.
#==============================================================================
        assert (temp > 70.0) and (temp < 3000.0), 'Gas temperature must be greater than 70K and less than 3000K.'
        if (wavemin == None):
            wavemin = amin(10000.0 / data['linecenter']) - 0.1
        if (wavemax == None):
            wavemax = amax(10000.0 / data['linecenter']) + 0.1
        ## First step: remove any data points that do not correspond to the primary isotope. (If we want to use isotopes,
        ## then we need to figure out mixing ratios.) For most HITRAN gases, the primary isotope is about 99% of the total
        ## atmospheric composition.
        okay = (data['I'] == 1)
        linecenters = array(data['linecenter'][okay])       ## line centers in wavenumbers
        linestrengths = array(data['S'][okay])*float(concentration.text())/100
        linewidths = array(data['gamma-air'][okay])* (1-float(concentration.text())/100) + array(data['gamma-self'][okay])*float(concentration.text())/100
        N_tempexps = array(data['N'][okay])     ## the temperature-dependent exponent for air-broadened linewidths
        nlines = alen(linecenters)
        Qratio = 1.0     # this is a placeholder for the ratio of total partition sums
        Epps = array(data['Epp'][okay])         ## the lower-energy-level energy (in cm^{-1})
        deltas = array(data['delta'][okay])     ## the "air pressure shift" (in cm^{-1} / atm)
        ## Convert the wavelengths (um) to wavenumbers (cm^{-1}). Create a spectrum linearly sampled in wavenumber (and
        ## thus nonuniformly sampled in wavelength).
        wavenumbers = linspace(10000.0/wavemax, 10000.0/wavemin, npts)
        waves = 10000.0 / wavenumbers
        xsec = zeros_like(wavenumbers)
        ## Define the list of channel boundary wavelengths.
        dk = wavenumbers[1] - wavenumbers[0]
        for i in arange(nlines):
            linecenter = linecenters[i]
            linestrength = linestrengths[i]
            linewidth = linewidths[i]
            N_tempexp = N_tempexps[i]
            Epp = Epps[i]
            delta = deltas[i]
            ## If the spectral line is well outside our region of interest, then ignore it.
            if (linecenter < amin(wavenumbers-0.5)):
                continue
            elif (linecenter > amax(wavenumbers+0.5)):
                continue
            ## If using a different temperature and pressure than the HITRAN default (296K and 1atm), then scale the
            ## linewidth by the temperature and pressure, adjust the linecenter due to pressure, and scale the
            ## linestrength.
            linecenter += delta * (pressure - 1.0) / pressure
            linewidth *= (pressure / 1.0) * pow(296.0/temp, N_tempexp)
            linestrength *= Qratio * exp(1.43877 * Epp * ((1.0/296.0) - (1.0/temp)))
            ## Note: the quantity sum(L * dk) should sum to "S"!
            L = lorentzian_profile(wavenumbers, linestrength, linewidth, linecenter)
            xsec += L
        if units.endswith('/mole'):
            xsec = xsec * 6.022E23
        elif units.endswith('.ppm'):
            xsec = xsec * 2.686E19
        if units.startswith('cm^2'):
            pass
        elif units.startswith('m^2'):
            xsec = xsec / 10000.0
        return(waves, xsec)

def lorentzian_profile(kappa, S, gamma, kappa0):
    '''
    Calculate a Lorentzian absorption profile.
    Parameters
    ----------
    kappa : ndarray
        The array of wavenumbers at which to sample the profile.
    S : float
        The absorption line "strength" (the integral of the entire line is equal to S).
    gamma : float
        The linewidth parameter of the profile.
    kappa0 : float
        The center position of the profile (in wavenumbers).
    Returns
    -------
    L : ndarray
        The sampled absorption profile.
    '''
    L = (S / pi) * gamma / ((kappa - kappa0)**2 + gamma**2)
    return(L)
    
if __name__ == "__main__": 
    app = QtGui.QApplication(sys.argv) 
    unit_um = True
    canvas = DataHolder()
    tabdialog = TabDialog() 
    tabdialog.show() 
    sys.exit(app.exec_()) 
    