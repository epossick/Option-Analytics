from wallstreet import Stock, Call, Put
import yfinance as yf
import statsmodels.api as sm
import scipy.stats as si
import pandas as pd
import math as m
import numpy as np




"""User Variables"""
# stock = ticker(s)
# period = time period
# interval = Sequence of time data
# valid period/interval variables: days = d, weeks = wk, months = mo, years = y, year-to-date = ytd, maximum = max 


"""Methods"""
# annual_stats()
# hedge_ratio()
# covariance_matrix
# correlation_matrix
#

"""Inputs"""
# Stats
stocks = "SPY QQQ"
period = "1y"
interval = "1d"


class Analytics():

    def __init__(self,ticker,period,interval):
        #Initialize Data
        self.split_tickers = ticker.split()
        split_tickers = ticker.split()
        self.period = period
        self.interval = interval
        if " " in ticker and len(list(set(split_tickers)))>1:
            tick=" ".join(list(set(split_tickers)))
            self.ticker = yf.Tickers(tick)
            self.history = self.ticker.history(period=self.period,interval=self.interval)
            self.clean = self.history.drop(["Open","High","Low","Volume","Dividends","Stock Splits"],axis=1)
            self.clean.columns = self.clean.columns.get_level_values(1)
        else:
            tick=" ".join(list(set(split_tickers)))
            self.ticker = yf.Ticker(tick)
            self.history = self.ticker.history(period=self.period,interval=self.interval).rename(columns={"Close": tick})
            self.clean = self.history.drop(["Open","High","Low","Volume","Dividends","Stock Splits"],axis=1)
        self.pct_changes = self.clean.pct_change().iloc[1:, :]
        
        #Initialize Dataframes
        periodmean = self.pct_changes.mean()
        periodvariance = self.pct_changes.var()
        periodstandard_dev = self.pct_changes.std()
        self.covariance_matrix = self.pct_changes.cov()
        self.period_mean = periodmean.to_frame().reset_index().rename(columns = {"index": "Ticker", 0: "Mean"})
        self.period_variance = periodvariance.to_frame().reset_index().rename(columns = {"index": "Ticker", 0: "Variance"})
        self.period_standard_dev = periodstandard_dev.to_frame().reset_index().rename(columns = {"index": "Ticker", 0: "Std"})
        self.covariance_matrix = self.pct_changes.cov()
        self.correlation_matrix = self.pct_changes.corr()
    
    def annual_stats(self):
        annual_stats=pd.DataFrame(columns=["Mean","Std"],index=self.period_mean["Ticker"])
        means=[]
        std = []
        if self.interval == "1d":
            time = 252
        elif self.interval == "1wk":
            time = 52
        elif self.interval == "1m":
            time = 12
        for i in self.period_mean["Mean"]:
            means.append((1 + i) ** time - 1)
        for i in self.period_standard_dev["Std"]:
            std.append((i * m.sqrt(time)))
        annual_stats["Mean"] = means
        annual_stats["Std"] = std
        return annual_stats

    def hedge_ratio(self):
        baseline = input("Which is the population: ")
        if self.split_tickers[0] == baseline:
            dependent = self.split_tickers[1]
        elif self.split_tickers[1] == baseline:
            dependent = self.split_tickers[0]
        correlation=self.correlation_matrix.iat[0,1]
        standard_deviations = self.annual_stats()["Std"]
        regression_coefficient = correlation * (standard_deviations[dependent] / standard_deviations[baseline])
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
        self.data = pd.DataFrame(rows,columns=["Ticker", "Strike", "Price","IV","Delta","Gamma","Theta","Vega","Flag", "Quantity"])
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
            dollar_delta = (data.at[i,"Delta"] * 100) * stock
            dollar_gamma = (data.at[i,"Gamma"] * 100) * (stock**2) / (100)
            levels1 = [self.split_tickers[i],levels_position[i]]
            levels.append(tuple(levels1))   
            labels.append([dollar_delta, dollar_gamma, data.at[i,"Theta"]*100, data.at[i,"Vega"]*100])           
        multi = pd.MultiIndex.from_tuples(levels)
        df = pd.DataFrame(labels, columns=cols, index=multi)
        #sum = df.sum().reset_index().rename(columns={0: })
        # Filtering Total data
        """for i in range(len(sum)):
            total.append(sum.at[i,0])"""
        df.loc['Net',:] = df.sum().values
        return df



s=Analytics(stocks,period,interval)
print(s.dollar_greeks())

