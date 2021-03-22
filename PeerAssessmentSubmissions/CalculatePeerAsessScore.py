import csv

import numpy as np
import glob
import pandas as pd

def read_file(fileName):
    """
    Read the input csv file and output the labels and the x data in np arrays
    :param fileName:
    :return: ys: N(N,); xs: M x N
    """
    with open(fileName) as f:
        reader = csv.reader(f, delimiter =',')
        data = list(reader)
        data = np.array(data)
    # print(data[0:2,:])
    # rows = data.shape[0]
    # cols = data.shape[1]
    # print('\n\n',fileName,'shape',data.shape, rows, cols)
    return data

def removeExtraStrings(data, startingIdx):
    for i in range(1,data.shape[0]): #skip row 0: header
        for j in range(startingIdx,startingIdx + 4): #for the rating columns; if
            # there
            # are other texts,
            # take the first digit only, ignore text.
            data[i, j] = data[i,j][0]
    # print('\ndone cleaning:\n ',data[0:2, :])
    return data

def getRelevantColumns(headerRow, cols, sorted=True):
    #usually in certain order, but one file is messed up, otherwise can
    # hardcode.
    if sorted:
        colIndex = [0,6,8,10,12,14,16] #select out the relevant columns: name,
        # 4 categories and comment
        return colIndex
    else:
        colStringIds = ['Name of your partner', 'statistical design', 'code',
                        'presentation','Communication', 'comments']
        colIndex = [0]
        for colId in colStringIds:
            for c in range(cols):
                if colId in headerRow[c]:
                    colIndex.append(c)
                    if len(colIndex) == 7:
                        return colIndex

def generateDataByName(data, startingIdx):
    dataByName = {} #initialize the dictionary
    for i in range(data.shape[0]):
        key = data[i,1].lower() #use last name as the key
        if '<' in key: #contains email address, remove
            key = key[:key.find('<')]
        if ',' in key: #last name first format
            key = key.split(',')[0]
        else: #last name last format
            key = key.split()[-1]
        key = key.strip() #remove extra white space
        key = key.replace(".", "") #remove the . symbol
        dataByName.setdefault(key, []).append(data[i,
                                              startingIdx:startingIdx+4].astype(int))
    return dataByName


def write_to_csv(totalPoints, outputFileName):
    with open(outputFileName, 'w') as f:
        for l in totalPoints:
            f.write(l[0] +',')
            f.write(str(l[1]) +'\n')


def processMidTermScore(fileList, offset):
    fulldata = [] #save all of them incase future processing is needed
    submitters = {}
    for fileIdx in range(len(fileList)):
        data = read_file(fileList[fileIdx])
        colIdx = getRelevantColumns(data[0, :], data.shape[1], fileIdx != 1)

        for data_row in range(1, data.shape[0]):
            lastnameKey = data[data_row, 0].split()[-1].lower()
            submitters.setdefault(lastnameKey, []).append('lab' + str(fileIdx + offset))

        data = data[:, colIdx]
        if fileIdx == 0:
            data = removeExtraStrings(data, 2)
        fileId = ['lab'+str(fileIdx+offset) for x in range(data.shape[0])]
        data = np.append(data, np.reshape(fileId, (len(fileId),1)), axis=1)
        if len(fulldata) == 0:
            fulldata = data
        else: #all following append, ignore the header row
            fulldata = np.append(fulldata, data[1:,:], axis=0)
    return fulldata, submitters


if __name__ == '__main__':
    fileList = glob.glob('Peer Assessment - Lab 0[2-4] Quiz Student Analysis '
                         'Report.csv')
    totalFiles = 3
    paneltyPts = 5
    fullPts = 20
    print(fileList)
    fulldata, submitters = processMidTermScore(fileList, 2)
    dataByName = generateDataByName(fulldata[1:,:], startingIdx=2)
    # print(dataByName)

    totalPts = {}
    for key,value in dataByName.items():
        pts = np.sum(value)
        if len(value) < totalFiles:
            diff = totalFiles - len(value)
            pts += diff * fullPts#partner has missing entry, default full points
        submissions = submitters[key]
        panelty = paneltyPts * (totalFiles - len(submissions))
        if (panelty > 0):
            print('panelty')
        totalPts[key] = np.sum(value) + (totalFiles - len(value)) * 20 - panelty
    # print(totalPts)
    totalPts = sorted(totalPts.items()) #sort by last name

    write_to_csv(totalPts, 'midterm_score_1.csv')

