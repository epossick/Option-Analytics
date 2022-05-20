from wallstreet import Stock, Call, Put
import yfinance as yf
import statsmodels.api as sm
import scipy.stats as si
import pandas as pd
import math as m
import numpy as np


"""Initial Inputs"""
stocks = "SPY"
period = "ytd"
interval = "1d"
key ="stats"
rounded = 3

"""Input Definitions"""
# stock = tickers (# of ticker inputs = # of positions  ex. "SPY SPY QQQ" for book with 2 SPY legs and 1 QQQ leg )
# period = time period
# interval = sequence of time data
    # valid period/interval: days = d, weeks = wk, months = mo, years = y, year-to-date = ytd, maximum = max 
# method = desired method within the class (include spaces for multiple methods ex. "stats cov" for annual_stats() and covariance_matrix)
# rounded = # of decimal places


"""Key Words for Methods"""

"""key"""        """Method"""
#  "stats"       annual_stats()
#  "hedge"       hedge_ratio()
#  "cov"         covariance_matrix
#  "correl"      correlation_matrix
#  "optdata"     option_data()
#  "$greek"      dollar_greeks()



class Analytics():

    def __init__(self,ticker,period,interval):
        # Initialize Data
        self.split_tickers = ticker.split()
        split_tickers = ticker.split()
        self.period = period
        self.interval = interval
        if " " in ticker and len(list(set(split_tickers)))>1:
            tick = " ".join(list(set(split_tickers)))
            self.ticker = yf.Tickers(tick)
            self.history = self.ticker.history(period=self.period,interval=self.interval)
            self.clean = self.history.drop(["Open","High","Low","Volume","Dividends","Stock Splits"],axis=1)
            self.clean.columns = self.clean.columns.get_level_values(1)
        else:
            tick = " ".join(list(set(split_tickers)))
            self.ticker = yf.Ticker(tick)
            self.history = self.ticker.history(period=self.period,interval=self.interval).rename(columns={"Close": tick})
            self.clean = self.history.drop(["Open","High","Low","Volume","Dividends","Stock Splits"],axis=1)
        self.pct_changes = self.clean.pct_change().iloc[1:, :]
        
        # Initialize Dataframes
        periodmean = self.pct_changes.mean()
        periodvariance = self.pct_changes.var()
        periodstandard_dev = self.pct_changes.std()
        self.period_mean = periodmean.to_frame().reset_index().rename(columns = {"index": "Ticker", 0: "E(r)"})
        self.period_variance = periodvariance.to_frame().reset_index().rename(columns = {"index": "Ticker", 0: "Real Var"})
        self.period_standard_dev = periodstandard_dev.to_frame().reset_index().rename(columns = {"index": "Ticker", 0: "Real Vol"})
        self.covariance_matrix = round(self.pct_changes.cov(), rounded)
        self.correlation_matrix = round(self.pct_changes.corr(), rounded)
    
    def annual_stats(self):
        # Convert period based return statistics to annual statistics
        annual_stats=pd.DataFrame(columns=["E(r)","Real Vol"],index=self.period_mean["Ticker"])
        means=[]
        std = []
        if self.interval == "1d":
            time = 252
        elif self.interval == "1wk":
            time = 52
        elif self.interval == "1m":
            time = 12
        for i in self.period_mean["E(r)"]:
            average = round((1 + i) ** time - 1, rounded)
            means.append(average)
        for i in self.period_standard_dev["Real Vol"]:
            stand_dev = round((i * m.sqrt(time)), rounded)
            std.append(stand_dev)
        annual_stats["E(r)"] = means
        annual_stats["Real Vol"] = std
        return annual_stats

    def hedge_ratio(self):
        # Approximate required ratio 
        baseline = input("What are you hedging?: ")
        if self.split_tickers[0] == baseline:
            hedge = self.split_tickers[1]
        elif self.split_tickers[1] == baseline:
            hedge = self.split_tickers[0]
        correlation=self.correlation_matrix.iat[0 , 1]
        standard_deviations = self.annual_stats()["Real Vol"]
        regression_coefficient = round(correlation * (standard_deviations[baseline] / standard_deviations[hedge]), rounded)
        return regression_coefficient
    
    def option_data(self):
        rows=[]
        for i in (self.split_tickers):
            self.flag = input("Flag: ")
            self.strike = float(input("strike: "))
            date = input("Date (m/d/yyyy): ") 
            self.quantity = int(input("Qty (minus for short): "))
            self.date = date.split("/")
            self.month = int(self.date[0])
            self.day = int(self.date[1])
            self.year = int(self.date[2])
            if self.flag == "call":
                self.opt = Call(i,self.day,self.month,self.year, self.strike)
            elif self.flag == "put":
                self.opt = Put(i,self.day,self.month,self.year, self.strike)
            self.option = self.opt.price
            rows.append([i, self.strike, self.option, self.opt.implied_volatility(), self.opt.delta()*self.quantity, \
                self.opt.gamma()*self.quantity, -1*self.opt.theta()*self.quantity, self.opt.vega()*self.quantity, self.flag, self.quantity])
        self.data = pd.DataFrame(rows,columns=["Ticker", "Strike", "Price", "IV", "Delta", "Gamma", "Theta", "Vega", "Flag", "Quantity"])
        return self.data

    def dollar_greeks(self):
        data = self.option_data()
        levels1=[]
        levels = []
        labels = []
        # Data in Dataframes
        levels_p = data.loc[:,"Quantity"].reset_index()
        levels_p["Quantity"] = np.where(levels_p["Quantity"] > 0, "Long", "Short")
        # Same Data in Lists 
        levels_position=levels_p["Quantity"].values.tolist()
        cols =("Delta", "Gamma", "Theta", "Vega")
        # Filter each position data
        for i in range(len(self.split_tickers)):
            stock = Stock(self.split_tickers[i]).price
            dollar_delta = (data.at[i , "Delta"] * 100) * stock
            dollar_gamma = (data.at[i , "Gamma"] * 100) * (stock ** 2) / (100)
            levels1 = [self.split_tickers[i],levels_position[i]]
            levels.append(tuple(levels1))   
            labels.append([round(dollar_delta, rounded), round(dollar_gamma, rounded), round(data.at[i , "Theta"] * 100, rounded), round(data.at[i , "Vega"] * 100, rounded)])           
        multi = pd.MultiIndex.from_tuples(levels)
        df = pd.DataFrame(labels, columns=cols, index=multi)
        df.loc['Net',:] = df.sum().values
        return df


# Instance for class
s=Analytics(stocks,period,interval)


# Assigns key to chosen method
if "stats" in key:
    print(s.annual_stats())
if "hedge" in key:
    print(s.hedge_ratio())
if "cov" in key:
    print(s.covariance_matrix)
if "correl" in key:
    print(s.correlation_matrix)
if "optdata" in key:
    print(s.option_data())
if "$greek" in key:
    print(s.dollar_greeks())



