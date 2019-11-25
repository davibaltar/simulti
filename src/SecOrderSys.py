
import wx
import numpy as np
import math

class TextOutput:
    def __init__ (self, zeta, omegaN, omegaD, tp, tr, ts, mp, sigma, theta, stateGraph, p1, p2):
        self.zeta = zeta
        self.omegaN = omegaN
        self.omegaD = omegaD
        self.tp = tp
        self.tr = tr
        self.ts = ts
        self.mp = mp
        self.sigma = sigma
        self.theta = theta
        self.stateGraph = stateGraph
        self.p1 = p1
        self.p2 = p2    

def compute(self, zeta, omegaN, amp):

    xGraph = 0
    yGraph = 0

    # Freq. natural amortecida
    omegaD = 0
    if zeta > 1:
        omegaD = omegaN * np.sqrt((zeta * zeta) - 1)
    else:
        omegaD = omegaN * np.sqrt(1 - (zeta * zeta))

    # Sigma
    sigma = zeta * omegaN

    # Angulo
    theta = 0
    if zeta != 0 and sigma != 0:
        theta = math.atan(omegaD / sigma)

    # Tempo de pico (Tp) para n = 1
    n = 1
    tp = 0
    if zeta < 1 and omegaD != 0:
        tp = (n * np.pi) / omegaD

    # Tempo de subida (Tr)
    tr = 0
    if zeta < 1 and omegaD != 0:
        tr = (np.pi - theta) / (omegaD)

    # Tempo de assentamento (Ts) 2%
    ts = 0
    if zeta != 0 and sigma != 0:
        ts = 3.912023 / sigma

    # Max. de sobressinal ou ultrapassagem (Conferir)
    mp = 0
    if zeta < 1:
        mp = amp * np.exp((-zeta * np.pi) / np.sqrt(1 - (zeta * zeta)))

    limDirGraph = 30

    t = np.arange(0.001, limDirGraph, 0.01)

    # Grafico de entrada R(s): Rampa, Degrau ou Impulso

    # Linha vermelha: Rampa (Implementar)
    # TODO

    # Linha vermelha: Degrau
    linha = amp * (t / t)

    # Linha vermelha: Impulso
    # TODO

    # Grafico Sinal Principal
    if zeta == 0:
        equacao = amp - (amp * np.cos(omegaN * t))
        xPoles = [0, 0]
        yPoles = [omegaN,-omegaN]
        p1Real = 0
        p1Imag = omegaN * np.sqrt(-((zeta * zeta) - 1))
        p2Real = 0
        p2Imag = -omegaN * (np.sqrt(-((zeta * zeta) - 1)))
        stateGraph = "Undamped" # "sem amortecimento"
    elif zeta > 0 and zeta < 1:
        equacao = amp - (amp * np.exp(-zeta * omegaN * t) * (np.cos(omegaN*t) + (zeta / (np.sqrt(1 - (zeta * zeta)))) * np.sin(omegaN * t)))
        xPoles = [-zeta * omegaN,-zeta * omegaN]
        yPoles = [omegaN * np.sqrt(-((zeta * zeta) - 1)), -omegaN * (np.sqrt(-((zeta * zeta) - 1)))]
        p1Real = -zeta*omegaN
        p1Imag = omegaN * np.sqrt(-((zeta * zeta) - 1))
        p2Real = -zeta * omegaN
        p2Imag = -omegaN * (np.sqrt(-((zeta * zeta) - 1)))
        stateGraph = "Underdamped" #"Sub-amortecido"
    elif zeta == 1:
        equacao = amp - (amp * np.exp((-omegaN * t) * (1 + omegaN * t)))
        xPoles = [-zeta * omegaN, -zeta * omegaN]
        yPoles = [0, 0]
        p1Real = -zeta * omegaN
        p1Imag = 0
        p2Real = -zeta * omegaN
        p2Imag = 0
        stateGraph = "Critically damped" # "Criticamente amortecido"
    elif zeta > 1:
        sg1 = (zeta + np.sqrt((zeta * zeta) - 1)) * omegaN
        sg2 = (zeta - np.sqrt((zeta * zeta) - 1)) * omegaN
        equacao = amp + (amp * (omegaN / (2 * np.sqrt(zeta*zeta - 1))) * ((np.exp(-sg1 * t) / sg1) - ( np.exp(-sg2 * t) / sg2)))
        xPoles = [(-zeta * omegaN) + np.sqrt(((zeta * zeta) - 1)), (-zeta * omegaN) - np.sqrt(((zeta * zeta) - 1))]
        yPoles = [0, 0]
        p1Real = (-zeta * omegaN) + np.sqrt(((zeta * zeta) - 1))
        p1Imag = 0
        p2Real = (-zeta * omegaN) - np.sqrt(((zeta * zeta) - 1))
        p2Imag = 0
        stateGraph = "Overdamped" # "Superamortecido"

    sign_im_p1 = ""
    sign_im_p2 = ""

    if p1Imag >= 0:
        sign_im_p1 = "+"
    elif p1Imag < 0:
        sign_im_p1 = " "
    
    if p2Imag >= 0:
        sign_im_p2 = "+"
    elif p2Imag < 0:
        sign_im_p2 = " "

    pole1 = str("{0:.2f}".format(p1Real))+" "+sign_im_p1+str("{0:.2f}".format(p1Imag))+"i"
    pole2 = str("{0:.2f}".format(p2Real))+" "+sign_im_p2+str("{0:.2f}".format(p2Imag))+"i"

    outputs = TextOutput(zeta, omegaN, omegaD, tp, tr, ts, mp, sigma, theta, stateGraph, pole1, pole2)
    return t, equacao, linha, xPoles, yPoles, outputs

def plotSecOrderSystem(self):
    amp = 0.0
    if self.radioImpulse.GetValue():
        amp = self.spinCtrlImpulse.GetValue()
    elif self.radioStep.GetValue():
        amp = self.spinCtrlStep.GetValue()
    else:
        amp = self.spinCtrlRamp.GetValue()

    xGraph, yGraph, linha, xPoles, yPoles, outputs = compute(self, float(self.spinCtrlZeta.GetValue()), float(self.spinCtrlOmega.GetValue()), amp) 
    
    # Main Graph
    self.plottedGraphs[0].set(xdata=xGraph, ydata=yGraph)
    self.plottedGraphs[1].set(xdata=xGraph, ydata=linha)
    self.mainGraph.set_xlim((0, 30))
    self.mainGraph.set_ylim((0, 2 * amp))

    # Poles Graph
    self.plottedGraphs[2].set(xdata=xPoles, ydata=yPoles)   
    self.polesGraph.set_xlim((-float(self.spinCtrlOmega.GetValue() + 1), 0))
    self.polesGraph.set_ylim((-(float(self.spinCtrlOmega.GetValue()) + 1), (float(self.spinCtrlOmega.GetValue()) + 1)))

    # Text Output
    self.outputParametersZeta.set_text(str("{0:.4f}".format(outputs.zeta)))
    self.outputParametersOmegaN.set_text(str("{0:.4f}".format(outputs.omegaN)))
    self.outputParametersOmegaD.set_text(str("{0:.4f}".format(outputs.omegaD)))
    self.outputParametersTp.set_text(str("{0:.4f}".format(outputs.tp)))
    self.outputParametersTr.set_text(str("{0:.4f}".format(outputs.tr)))
    self.outputParametersTs.set_text(str("{0:.4f}".format(outputs.ts)))
    self.outputParametersMp.set_text(str("{0:.4f}".format(outputs.mp)))
    self.outputParametersSigma.set_text(str("{0:.4f}".format(outputs.sigma)))
    self.outputParametersTheta.set_text(str("{0:.4f}".format(outputs.theta)))
    self.outputParametersState.set_text(outputs.stateGraph)
    self.outputParametersP1.set_text("Pole 1: " + outputs.p1)
    self.outputParametersP2.set_text("Pole 2: " + outputs.p2)

    return 0
