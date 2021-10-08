#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  6 17:56:23 2021

@author: alok

"""

def extract_symmetry(emdid):
    from selenium import webdriver
    from selenium.webdriver.firefox.options import Options
    from bs4 import BeautifulSoup as soup
    import pandas as pd
    
    ## Open browser in headless 
    opts = Options()
    opts.add_argument("--headless")
    driver = webdriver.Firefox(options=opts)
    
    ## Script to web scrap and get symmetry information
    
    emdb_url = "https://www.emdataresource.org/EMD-"+str(emdid)
    
    driver.get(emdb_url)
    page_html = driver.page_source
    driver.close()
    
    ## Parse using pandas
    
    dataframes = pd.read_html(page_html)
    
    symmetry = "unknown"

    for df in dataframes:
        if 0 in list(df.keys()):
            if "Imposed Symmetry" in df[0].values:
                col = df[0]
                symmetry = df[col=="Imposed Symmetry"][1].values[0]
                print("Found symmetry is: "+symmetry)
    
    return symmetry
