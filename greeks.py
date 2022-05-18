from subprocess import call
from requests import put
from wallstreet import Stock
from math import sqrt, log, exp
import scipy.stats as stat
from scipy.stats import norm

"""This program Provides Greeks for """
class Greeks:
    #ticker = stock ticker
    #flag = put/call
    #strike = option strike
    #time = time to maturity (days of CALENDAR year)
    #iv = implied volatility (annual %)
    #rf = risk - free rate (annual %)
    #div = dividend yield (annual %)
    def __init__(self, flag, ticker, strike, time, iv, rf, div):
        self.flag = flag
        self.stock = float(Stock(ticker).price)
        self.strike = float(strike)
        self.days = 365
        self.T = float(time) / self.days
        self.sigma = float(iv)
        self.rf = float(rf)
        self.div = float(div)
        self.sigmaT = self.sigma * sqrt(self.T)
        self.d1 = (log(self.stock / self.strike) + (self.rf - self.div+(self.sigma**2/2))*self.T) / self.sigmaT
        self.d2 = self.d1 - self.sigmaT
        self.option = self.Option()
        # First Order Greeks
        self.delta = round(self.Delta(),5)
        self.theta = round(self.Theta(),5)
        self.vega = round(self.Vega(),5)
        self.rho = round(self.Rho(),5)  
        self.phi = round(self.Phi(),5)
        self.Lambda = self.Elasticity()
        # Second Order Greeks
        self.gamma = round(self.Gamma(),5)
        self.vanna = round(self.DdeltaDvol(),5)
        self.charm = round(self.DdeltaDtime(),5)
        self.volga = round(self.DvegaDvol(),5)
        self.veta  = round(self.DvegaDtime(),5)
        # Third Order Greeks
        self.speed = round(self.DgammaDspot(),5)
        self.zomma = round(self.DgammaDvol(),5)
        self.color = round(self.DgammaDtime(),5)
        self.ultima = round(self.DvommaDvol(),5)
        # Dollar Greeks
        self.Ddelta = round(self.Dollar_Delta(),5)
        self.Dgamma = round(self.Dollar_Gamma(),5)

    """Option price"""
    def Option(self):
        dfq = exp(-self.div * self.T)
        df = exp(-self.rf * self.T)
        if self.flag == "c":
            price = self.stock * dfq * norm.cdf(self.d1) - self.strike * df * norm.cdf(self.d2)
        if self.flag =="p":
            price = self.strike * df * norm.cdf(-self.d2) - self.stock * dfq * norm.cdf(-self.d1)
        return price

    ####################################################################
    ######################### 1st Order Greeks #########################
    ####################################################################

    def Delta(self):
        dfq = exp(-self.div * self.T)
        if self.flag == "c":
            return dfq * norm.cdf(self.d1) 
        if self.flag == "p":
            return dfq * (norm.cdf(self.d1) - 1) 
    
    #Theta in calendar days
    def Theta(self):
        dfq = exp(-self.div * self.T)
        df = exp(-self.rf * self.T)
        if self.flag == "c":
            return (1 / self.days) * ((-1 * (self.stock * self.sigma * dfq) / (2 * sqrt(self.T)) * norm.pdf(self.d1)) \
                - self.rf * self.strike * df * norm.cdf(self.d2) + self.div * self.stock * dfq * norm.cdf(self.d1))
        if self.flag == "p":
            return (1 / self.days) * ((-1 * (self.stock * self.sigma * dfq) / (2 * sqrt(self.T)) * norm.pdf(-self.d1)) \
                + self.rf * self.strike * df * norm.cdf(-self.d2) - self.div * self.stock * dfq * norm.cdf(-self.d1))
        
    def Vega(self):
        dfq = exp(-self.div * self.T)
        return (1 / 100) * (self.stock * dfq * sqrt(self.T) * norm.pdf(self.d1))
    
    def Rho(self):
        df = exp(-self.rf * self.T)
        if self.flag == "c":
            return (1 / 100) * (self.strike * self.T * df * norm.cdf(self.d2))
        if self.flag == "p":
            return (-1 / 100) * (self.strike * self.T * df * norm.cdf(-self.d2))
    
    def Phi(self):
        dfq = exp(-self.div * self.T)
        if self.flag == "c":
            return -(1 / 100) * (self.stock * self.T * dfq * norm.cdf(self.d1))
        if self.flag == "p":
            return -(1 / 100) * (self.stock * self.T * dfq * norm.cdf(-self.d1))

    def Elasticity(self):
        return self.delta * (self.stock / self.option)

    ####################################################################
    ######################### 2nd Order Greeks #########################
    ####################################################################

    def Gamma(self):
        dfq = exp(-self.div * self.T)
        return dfq / (self.stock * self.sigmaT) * norm.pdf(self.d1)

    def DdeltaDvol(self):
        dfq = exp(-self.div * self.T)
        if self.flag == "c":
            return (1 / 100) * (((self.d2 * dfq * norm.pdf(self.d1)) * -1) / (self.sigma))
        if self.flag == "p":
            return (1 / 100) * (((self.d2 * dfq * norm.pdf(-self.d1)) * -1) / (self.sigma))

    def DdeltaDtime(self):
        dfq = exp(-self.div * self.T)
        if self.flag == "c":
            return (1 / 365) * (self.div * dfq * norm.cdf(self.d1) + dfq * norm.pdf(self.d1) \
                * ((self.d2 / (2 * self.T)) - ((self.rf - self.div) / (self.sigmaT))))
        if self.flag == "p":
            return (1 / 365) * (- 1 * self.div * dfq * norm.cdf(-self.d1) + dfq * norm.pdf(self.d1) \
                * ((self.d2 / (2 * self.T)) - ((self.rf - self.div) / (self.sigmaT))))
             
    def DvegaDvol(self):
        dfq = exp(-self.div * self.T)
        if self.flag == "c":
            return (1 / 100) * (self.d1 * self.d2 * self.stock * dfq * norm.pdf(self.d1) * self.T) / (self.sigma)
        if self.flag == "p":
            return (1 / 100) * (self.d1 * self.d2 * self.stock * dfq * norm.pdf(-self.d1) * self.T) / (self.sigma)

    def DvegaDtime(self):
        dfq = exp(-self.div * self.T)
        return (1 / (self.days * 100)) * (-1 * self.stock * dfq * norm.pdf(self.d1) * sqrt(self.T) \
            * (self.div + (((self.rf - self.div) * self.d1) / self.sigmaT) - ((1 + self.d1 * self.d2)/(2 * self.T))))

    ####################################################################
    ######################### 3rd Order Greeks #########################
    ####################################################################

    def DgammaDspot(self):
        dfq = exp(-self.div * self.T)
        return (-1 * dfq) * (norm.pdf(self.d1) / (self.stock ** 2 * self.sigmaT)) * (self.d1 / self.sigmaT + 1)
    
    def DgammaDvol(self):
        dfq = exp(-self.div * self.T)
        return (1 / 100) * (dfq * (norm.pdf(self.d1) * (self.d1 * self.d2 - 1)) / (self.stock * self.sigma ** 2 * sqrt(self.T)))
    
    def DgammaDtime(self):
        dfq = exp(-self.div * self.T)
        return (1 / self.days) * (-1 * dfq * norm.pdf(self.d1) / (2 * self.stock * self.T * self.sigmaT) \
             * (2 * self.div * self.T + 1 + (2 * (self.rf - self.div) * self.T - self.d2 * self.sigmaT / self.sigmaT) * self.d1))

    def DvommaDvol(self):
        return  (1 / 100) * (-self.vega / (self.sigma ** 2) * (self.d1 * self.d2 * (1 - self.d1 * self.d2) + self.d1 ** 2 + self.d2 ** 2))
        
    ####################################################################
    ######################### Dollar Greeks ############################
    ####################################################################

    def Dollar_Delta(self):
        return (self.delta * 100) * self.stock

    def Dollar_Gamma(self):
        return (self.gamma * 100) * (self.stock ** 2) / (100)



    