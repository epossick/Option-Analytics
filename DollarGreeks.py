from wallstreet import Stock, Call, Put
import numpy as np
import pandas as pd

class DollarGreeks:
    # ticker = stock ticker
    # flag = call/put
    # strike = option strike
    # date = maturity date

    def __init__(self, ticker, flag, strike, date):
        self.ticker = ticker
        self.stock=Stock(self.ticker).price
        self.flag = flag
        self.strike = float(strike)
        self.date = date.split("/")
        self.month = int(self.date[0])
        self.day = int(self.date[1])
        self.year = int(self.date[2])
        if self.flag == "call":
            option = Call(self.ticker,self.day,self.month,self.year, self.strike)
        elif self.flag == "put":
            option = Put(self.ticker,self.day,self.month,self.year, self.strike)
        self.option = option.price
        self.IV = option.implied_volatility()
        self.delta = option.delta()
        self.gamma = option.gamma()
        self.theta = option.theta()
        self.vega = option.vega()
    # Theta and Vega are already quoted in Dollars
    def dollar_greeks(self):
        greeks=pd.DataFrame()
        dollar_delta = (self.delta*100) * self.stock
        dollar_gamma = (self.gamma*100) * (self.stock**2) / (100)
        greeks['Ticker'] = [self.ticker]
        greeks["Dollar Delta"] = [dollar_delta]
        greeks["Dollar Gamma"] = [dollar_gamma]
        return greeks
"""ticker = input("What is the ticker: ")
flag = input("What is the flag (call/put): ")
date = input("What is the date (m/d/y): ")
strike = input("What is the strike: ")"""
ticker="TSLA"
flag="call"
date="3/11/2022"
strike="840"

option1=Greeks(ticker,flag,strike,date)
print(option1.dollar_greeks())
