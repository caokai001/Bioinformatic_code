# 两个文档链接
#1.http://jaspar.genereg.net/api/overview
#2.http://jaspar.genereg.net/api/clients

import coreapi
import requests
import os

# Initialize a client & load the schema document
client = coreapi.Client()
schema = client.get("http://jaspar.genereg.net/api/v1/docs")
 
# Interact with the API endpoint
action = ["matrix", "list"]

def modify_params(name='FOXA1',tax_id="9606"):
   
    '''
    主要是修改params 参数，用于api查询.主要参数如下：
    
    name  :   TF name
    tax_id:   specise
    
    '''
    params = {
    "collection": 'CORE',
    "name": name,
    "tax_group": 'vertebrates',
    "tax_id": tax_id,  #物种
    "version": 'latest',
    "release": '2018',
}
    return params

#def download_svg(url):
#    '''
#    下载motif log 图片
#    
#    '''
#    print(file_name)
#    #file_name=result["results"][0]['name']+"_"+result["results"][0]['url'].split("/")[-2]
#    #path=save_dir+"\\"+file_name
#
#    print(url)
#    try:
#        r=requests.get(url)
#        r.raise_for_status()
#        with open(file_name+".svg",'wb') as f:  # 二进制文件
#            f.write(r.content)
#    except:
#        raise "error exit !"



     
def modify_matrix_format(matrix_raw_url,fmt="jaspar"):
    '''
    修改matrix_raw_url ,添加format.默认是jaspar.支持很多格式，
    详见：http://jaspar.genereg.net/api/overview
    available data fromats are; json, jsonp, api, yaml, jaspar, transfac, pfm
    
    '''
    
    matrix_url=matrix_raw_url+"?format={}".format(fmt)
    return matrix_url       
    

    
#def download_matrix(url):
#    '''
#    下载特定格式的matrix,如jasper 
#    
#    '''
#
#    
#    try:
#        r=requests.get(url)
#        r.raise_for_status()
#        with open(file_name+"_matrix.txt",'w') as f:  # 二进制文件
#            f.write(r.text)
#        
#    except:
#        raise "error exit !"
    

    
    
 

def download_single(name,tax_id="9606"):
    '''
    download only one matrix and picture
    
    '''
    print(name)
    params=modify_params(name=name)
    print(params)
    result = client.action(schema, action, params=params)
    
    file_name=result["results"][0]['name']+"_"+result["results"][0]['url'].split("/")[-2]
    
    
    print(file_name)
    # 建立文件夹
    try:
        os.mkdir(file_name)
    except FileExistsError:
        pass
    
    os.chdir(file_name)
    
    # 下载seq_log svg
    seqlog_url=result["results"][0]["sequence_logo"]
    print(seqlog_url)
    
    try:
        print(seqlog_url)
        r=requests.get(seqlog_url)
        r.raise_for_status()
        print(r)
        with open(file_name+".svg",'wb') as f:  # 二进制文件
            print("write ok")
            f.write(r.content)
    except:
        raise "error exit seq!"
    
    #download_svg(seqlog_url)
    # 下载matrix
    matrix_raw_url=result["results"][0]["url"]
    matrix_url=modify_matrix_format(matrix_raw_url,fmt="jaspar")
    try:
        r=requests.get(matrix_url)
        r.raise_for_status()
        with open(file_name+"_matrix.txt",'w') as f:  # 二进制文件
            f.write(r.text)
        
    except:
        raise "error exit !"
    
#    print(matrix_url)
#    download_matrix(matrix_url)



    
def main(L):
    '''
    判断data 是文件还是列表list,tuple
    
    
    '''
#    L=[]
#    if os.path.isfile(data):
#        with open(data) as f:
#            L=[ line.strip() for line in f.readlines()]
#    elif isinstance(data,(list,tuple)):
#        L=data
#        
#    print(L)
    
    os.chdir("C:\\Users\\hp\\Desktop\\python并行")
    for i in L:
        download_single(name=i)
        os.chdir("C:\\Users\\hp\\Desktop\\python并行")

L=['CTCF', 'FOXA1', 'AR']
main(L)
