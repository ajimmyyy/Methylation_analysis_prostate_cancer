from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.common.action_chains import ActionChains
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.support.select import Select
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
import pandas as pd
import numpy as np
import time

gene_list_dir = "C:/Users/acer/Desktop/special_subject/code/gene.txt"
fn_dmp = "DMP_hyper_with_GO_all.csv"
fn_o = "DMP_hyper_with_GO_all.csv"
fn_o_func = "function_surface_MF.txt"
time_choose = 1000

def find_func(row, dict):
    if row["gene"] in dict.keys():
        return ",".join(str(x) for x in dict[row["gene"]]) + ","
    else:
        return None
    
options = webdriver.ChromeOptions()
options.add_experimental_option('detach', True)

service = Service(executable_path = "C:/Users/acer/Desktop/special_subject/code/download_go_term/chromedriver.exe")
chrome = webdriver.Chrome(service=service, options=options)

# 輸入基因
chrome.get("https://david.ncifcrf.gov/tools.jsp")
chrome.implicitly_wait(20)

# 主窗口
main_window_handle = chrome.current_window_handle
now_window_handle = main_window_handle

chrome.find_element(By.XPATH, '//*[@id="tablist"]/li[1]/a').click()
chrome.find_element(By.XPATH, '//*[@id="divUpload"]/table/tbody/tr[9]/td/p/input').send_keys("C:/Users/acer/Desktop/special_subject/code/gene.txt")
Select(chrome.find_element(By.XPATH, '//*[@id="Identifier"]')).select_by_value("OFFICIAL_GENE_SYMBOL")
chrome.implicitly_wait(3)
chrome.find_element(By.XPATH, '//*[@id="speciesSelect"]').send_keys("Homo sapiens")
ActionChains(chrome).move_to_element(chrome.find_element(By.XPATH, '//*[@id="speciesSelectDiv"]/span/div/div/div[1]')).click().perform()
chrome.find_element(By.XPATH, '//*[@id="divUpload"]/table/tbody/tr[17]/td/table/tbody/tr[1]/td[2]/input').click()
chrome.find_element(By.XPATH, '//*[@id="divUpload"]/table/tbody/tr[18]/td/input').click()

chrome.find_element(By.XPATH, '/html/body/table[2]/tbody/tr/td[3]/table/tbody/tr/td/table/tbody/tr/td/table/tbody/tr/td/table/tbody/tr[5]/td[2]/a').click()
chrome.find_element(By.XPATH, '//*[@id="summaryTree"]/li[3]/a').click()

# BP
# chrome.find_element(By.XPATH, '//*[@id="summaryTree"]/li[3]/ul/li/table/tbody/tr[7]/td[4]/input').click()
# CC
# chrome.find_element(By.XPATH, '//*[@id="summaryTree"]/li[3]/ul/li/table/tbody/tr[15]/td[4]/input').click()
# MF
chrome.find_element(By.XPATH, '//*[@id="summaryTree"]/li[3]/ul/li/table/tbody/tr[23]/td[4]/input').click()

# 切換至副視窗
WebDriverWait(chrome, 10).until(lambda chrome: len(chrome.window_handles) > 1)
all_window_handles = chrome.window_handles

for handle in all_window_handles:
    if handle != main_window_handle:
        chrome.switch_to.window(handle)
        now_window_handle = chrome.current_window_handle
        break

chrome.find_element(By.XPATH, '//*[@id="row"]/thead/tr/th[6]/a').click()
img_elements = chrome.find_elements(By.XPATH, "//tr[@class='odd' or @class='even']//img")
func_elements = chrome.find_elements(By.XPATH, "//tr[@class='odd' or @class='even']/td[position()=3]")
func_list = [td.text for td in func_elements]
gene_list = []
for i in range(time_choose):
    if i >= len(img_elements):
        break
    img_element = img_elements[i]
    img_element.click()
    all_window_handles = chrome.window_handles
    for handle in all_window_handles:
        if handle != main_window_handle and handle != now_window_handle:
            chrome.switch_to.window(handle)
            break

    gene_names = chrome.find_elements(By.XPATH, "//tr[@class='odd' or @class='even']/td[position()=1]")
    gene_list_per_func = [td.text for td in gene_names]
    gene_list.append(gene_list_per_func)
    chrome.close()
    chrome.switch_to.window(now_window_handle)
chrome.close()

# 建立gene對功能的字典
gene_to_func = {}
for i in range(len(gene_list)):
    for j in gene_list[i]:
        if j not in gene_to_func:
            gene_to_func[j] = [i]
        else:
            gene_to_func[j].append(i)

# 對照表
with open(fn_o_func, "w") as file:
    for idx, item in enumerate(func_list):
        file.write(f"{idx}: {item}\n")

# csv檔
data_dmp_df = pd.read_csv(fn_dmp)
# data_dmp_df["BP"] = data_dmp_df.apply(find_func, dict = gene_to_func, axis = 1)
# data_dmp_df["CC"] = data_dmp_df.apply(find_func, dict = gene_to_func, axis = 1)
data_dmp_df["MF"] = data_dmp_df.apply(find_func, dict = gene_to_func, axis = 1)
data_dmp_df.to_csv(fn_o, sep=',', encoding='utf-8', index=False)

# 關閉
chrome.quit()