#include <iostream>
#include<stdlib.h>
#include <fstream>
#include <vector>
#include <queue>
#include <algorithm>
#include <math.h>
#include <string>
using namespace std;
typedef struct chipUnitBlock{
    int m1Blocked;
    int m2Blocked;
    int via;//0-no via 1-via 2- beside via block
    int cellNumber;
    int m1NetNumber;
    int m2NetNumber;
    int m1LeeNumber;
    int m2LeeNumber;
    int cellSurrounding;
    int cellTerminal;
    int m1NetSurrounding;
    int m2NetSurrounding;
    int x;
    int y;
}chipUnitBlock;
typedef struct coOrdinates{
    int x;
    int y;
}coOrdinates;
typedef struct emptyCellInfo{
    int x;
    int y;
    int distance;
}emptyCellInfo;

typedef struct cellcoOrdinates{
    coOrdinates referencePoint;
    coOrdinates t1;
    coOrdinates t2;
    coOrdinates t3;
    coOrdinates t4;
    int orientation;
}cellcoOrdinates;
vector<cellcoOrdinates> cellCoOrdinateVector;
typedef struct netInfo{
    int netNumber;
    int netPriority;
    int sourceCellNumber;
    int sourceCellTerminalNumber;
    int targetCellNumber;
    int targetCellTerminalNumber;
    int routingMetal;
}netInfo;
typedef struct cellWeightInfo{
    int cellNumber;
    int cellWeight;
    bool blocked;
}cellWeightInfo;
typedef struct cellInfo{
    int x;
    int y;
    int cellWeightVectorIndex;
    bool blocked;
    int cellNumber;
}cellInfo;
vector< vector<chipUnitBlock> > chipGridMatrix;
vector< vector<int> > cellConnectivityMatrix;
vector< vector<cellInfo> > cellPlacementMatrix;
vector<cellWeightInfo> cellWeightVector;
vector<cellInfo> cellPositionVector;
vector<netInfo> netTerminalVector;
queue<chipUnitBlock> unVisitedQueue;

void saveNetInformation(ifstream &benchMarkFile,int &noOfCells,int &noOfNets);
void convertPlacementMatrixToChipGridMatrix(int noOfCells,int cellToBeFlipped,int rowSpacing,int columnSpacing);
void printChipGridMatrix();
void printChipGridMatrixElement(int i,int j);
void printPlacementMatrix();
bool isCompletelyBlocked(chipUnitBlock tempChipUnitBlock);
bool isPointInGridBounds(int x, int y);
void blockCellsAroundTheBendForMetal1(int netNumber,chipUnitBlock bendBlock);
void unBlockCellsAroundTheBendForMetal1(int netNumber,chipUnitBlock bendBlock,chipUnitBlock source,chipUnitBlock target);
void blockCellsAroundTheBendForMetal2(int netNumber,chipUnitBlock bendBlock);
void unBlockCellsAroundTheBendForMetal2(int netNumber,chipUnitBlock bendBlock,chipUnitBlock source,chipUnitBlock target);
void insertViaAtGivenLocation(int netNumber,chipUnitBlock via);
bool checkIfViaCanBeInserted(chipUnitBlock via);
void deleteViaAtGivenLocation(chipUnitBlock viaBlock);
int generateLeeNumbersFromSourceToTerminal(chipUnitBlock source,chipUnitBlock target,int metalType);
int backTraceNetPath(int netNumber,chipUnitBlock source,chipUnitBlock target);
void clearLeeNumbersFromSourceToTerminal(int netNumber,chipUnitBlock source,chipUnitBlock target);
int performLeeAlgorithmForTheNets(int noOfNets,int &failedMethod,int metalType);
void placeCellsAtRandomLocations(int noOfCells);
coOrdinates findTargetEmptyLocation(coOrdinates source);
bool compareEmptyCellInfo(emptyCellInfo cell1,emptyCellInfo cell2);
int getWireLengthEstimate();
void runForceDirectedPlacementAlgorithm(int noOfCells);
bool compareCellInfo(cellWeightInfo node1,cellWeightInfo node2);
coOrdinates getTerminalCoOrdinates(cellcoOrdinates cell,int terminalNumber);
void convertChipGridMatrixToMagicFile(ofstream &outputMagicFile);
bool compareNetInfo(netInfo net1,netInfo net2);
bool sortByNetNumber(netInfo net1,netInfo net2);
void prioritizeNets();
float areaOfRect(coOrdinates rect1,coOrdinates rect2,coOrdinates rect3,coOrdinates rect4);
float areaOfTriangle(coOrdinates a,coOrdinates b,coOrdinates p);
bool checkIfAPointIsInsideARectangle(coOrdinates rect1,coOrdinates rect2,coOrdinates rect3,coOrdinates rect4,coOrdinates point);
void flipACell(int cellNumber,int orientation);
int getFinalOrientation(int orientation);
void clearNetSurroundingNumber(int metalType,int netNumber,chipUnitBlock currentBlock);
void updateNetSurroundingNumber(int metalType,int netNumber,chipUnitBlock currentBlock);
bool isBlockInSurrounding(chipUnitBlock currentBlock,chipUnitBlock source,chipUnitBlock target);

int main(int argumentCount, char* fileName[]){
    int noOfCells = 0;
    int noOfNets = 0;
    ifstream benchMarkFile(fileName[1]);
    ofstream outputMagicFile(fileName[2]);
    int unRoutableNet = -1;
    //int leeLimit = 50;
    int failedMethod = 0;
    int wireLength = 0;
    int rowSpacing = 12;
    int columnSpacing = 12;
    if(benchMarkFile.is_open()){
        saveNetInformation(benchMarkFile,noOfCells,noOfNets);//common
        //****** placement begins**********
        placeCellsAtRandomLocations(noOfCells);
        wireLength = getWireLengthEstimate();
        cout<<"Initial wireLength: "<< wireLength <<endl;
        runForceDirectedPlacementAlgorithm(noOfCells);
        //printPlacementMatrix();
        wireLength = getWireLengthEstimate();
        cout<<"Wire length after final placement: "<< wireLength <<endl;
        //****** placement ends**********
        if(noOfNets>500){ // for better routability
            rowSpacing = 30;
            columnSpacing = 30;
        }else if(noOfNets>600){
            rowSpacing = 40;
            columnSpacing = 40;
        }
        convertPlacementMatrixToChipGridMatrix(noOfCells,-1,rowSpacing,columnSpacing);//common
        //****** routing begins**********
        prioritizeNets();
        unRoutableNet = performLeeAlgorithmForTheNets(noOfNets,failedMethod,1);
        //****** routing ends**********
        convertChipGridMatrixToMagicFile(outputMagicFile);//common
    }else{
        cout<<"File could not be opened.\n"<<endl;
    }
    return 0;
}
void saveNetInformation(ifstream &benchMarkFile,int &noOfCells,int &noOfNets){
    benchMarkFile>>noOfCells;
    benchMarkFile>>noOfNets;
    cellConnectivityMatrix = vector< vector<int> >(noOfCells);
    int cellConnectivityMatrixSize = 0;
    int netNumber = 0;
    int sourceCellNumber =0;
    int sourceCellTerminalNumber =0;
    int targetCellNumber = 0;
    int targetCellTerminalNumber = 0;
    netInfo netInfoStruct;
    cellWeightInfo currentCellWeightInfo;
    while(benchMarkFile >> netNumber >> sourceCellNumber >> sourceCellTerminalNumber >> targetCellNumber >> targetCellTerminalNumber){
        netInfoStruct.netNumber = netNumber;
        netInfoStruct.netPriority = 0;
        netInfoStruct.sourceCellNumber = sourceCellNumber;
        netInfoStruct.sourceCellTerminalNumber = sourceCellTerminalNumber;
        netInfoStruct.targetCellNumber = targetCellNumber;
        netInfoStruct.targetCellTerminalNumber = targetCellTerminalNumber;
        netInfoStruct.routingMetal = -1;
        netTerminalVector.push_back(netInfoStruct);
        cellConnectivityMatrix[sourceCellNumber-1].push_back(targetCellNumber-1);
        cellConnectivityMatrix[targetCellNumber-1].push_back(sourceCellNumber-1);
    }
    cellConnectivityMatrixSize = cellConnectivityMatrix.size();
    for(int i=0;i<cellConnectivityMatrixSize;i++){
        currentCellWeightInfo.cellNumber = i;
        currentCellWeightInfo.cellWeight = cellConnectivityMatrix[i].size();
        currentCellWeightInfo.blocked = false;
        cellWeightVector.push_back(currentCellWeightInfo);
    }
}

void placeCellsAtRandomLocations(int noOfCells){
    int cellMatrixSize = ceil(sqrt(noOfCells));
    cellPlacementMatrix = vector< vector<cellInfo> >(cellMatrixSize,vector<cellInfo>(cellMatrixSize));
    cellPositionVector.resize(noOfCells);
    int cellNumber = 0;
    cellWeightInfo currentCellWeightInfo;
    cellInfo sourceCellInfo;
    for(int i=0;i<cellMatrixSize;i++){
        for(int j=0;j<cellMatrixSize;j++){
            if(cellNumber<noOfCells){
                sourceCellInfo.blocked = false;
                sourceCellInfo.cellWeightVectorIndex = -1;
                sourceCellInfo.x = i;
                sourceCellInfo.y = j;
                sourceCellInfo.cellNumber = cellNumber;
                cellPlacementMatrix[i][j] = sourceCellInfo;
                cellPositionVector[cellNumber] = sourceCellInfo;
                cellNumber++;
            }else{
                sourceCellInfo.x = i;
                sourceCellInfo.y = j;
                sourceCellInfo.cellNumber = -1;
                sourceCellInfo.cellWeightVectorIndex = -1;
                sourceCellInfo.blocked  = false;
                cellPlacementMatrix[i][j] = sourceCellInfo;
            }
        }
    }
    sort(cellWeightVector.begin(),cellWeightVector.end(),compareCellInfo);
    for(int i=0;i<noOfCells;i++){
        currentCellWeightInfo = cellWeightVector[i];
        int currentCellNumber = currentCellWeightInfo.cellNumber;
        sourceCellInfo = cellPositionVector[currentCellNumber];
        sourceCellInfo.cellWeightVectorIndex = i;
        cellPositionVector[currentCellNumber] = sourceCellInfo;
        cellPlacementMatrix[sourceCellInfo.x][sourceCellInfo.y] = sourceCellInfo;
    }
}

int getWireLengthEstimate(){
    int cellConnectivityMatrixSize = cellConnectivityMatrix.size();
    int currentCellConnections = 0;
    coOrdinates source;
    coOrdinates target;
    int wireLength = 0;
    for(int i=0;i<cellConnectivityMatrixSize;i++){
        currentCellConnections = cellConnectivityMatrix[i].size();
        source.x = cellPositionVector.at(i).x;
        source.y = cellPositionVector.at(i).y;
        for(int j=0;j<currentCellConnections;j++){
            target.x = cellPositionVector.at(j).x;
            target.y = cellPositionVector.at(j).y;
            wireLength += abs(source.x-target.x)+abs(source.y-target.y);
        }
    }
    return wireLength/2;
    /*int netTerminalVectorSize = netTerminalVector.size();
    netInfo currentNetInfo;
    coOrdinates source;
    coOrdinates target;
    int wireLength = 0;
    for(int i=0;i<netTerminalVectorSize;i++){
        currentNetInfo = netTerminalVector.at(i);
        source = getTerminalCoOrdinates(cellCoOrdinateVector.at(currentNetInfo.sourceCellNumber-1),currentNetInfo.sourceCellTerminalNumber);
        target = getTerminalCoOrdinates(cellCoOrdinateVector.at(currentNetInfo.targetCellNumber-1),currentNetInfo.targetCellTerminalNumber);
        wireLength += abs(source.x-target.x)+abs(source.y-target.y);
    }
    return wireLength;*/
}

coOrdinates findTargetEmptyLocation(coOrdinates source){
    int cellMatrixSize = cellPlacementMatrix.size();
    vector<emptyCellInfo> emptyCellVector;
    emptyCellInfo currentCellInfo;
    coOrdinates target;
    for(int i=0;i<cellMatrixSize;i++){
        for(int j=0;j<cellMatrixSize;j++){
            if(cellPlacementMatrix[i][j].cellNumber == -1){
                currentCellInfo.x = i;
                currentCellInfo.y = j;
                //currentCellInfo.distance = abs(source.x-currentCellInfo.x)+abs(source.y-currentCellInfo.y);
                currentCellInfo.distance = sqrt((source.x-currentCellInfo.x)*(source.x-currentCellInfo.x) +(source.y-currentCellInfo.y)*(source.y-currentCellInfo.y));
                emptyCellVector.push_back(currentCellInfo);
            }
        }
    }
    if(emptyCellVector.size()>0){
        sort(emptyCellVector.begin(),emptyCellVector.end(),compareEmptyCellInfo);
        target.x = emptyCellVector.at(0).x;
        target.y = emptyCellVector.at(0).y;
        return target;
    }
    return source;

}

bool compareEmptyCellInfo(emptyCellInfo cell1,emptyCellInfo cell2){
    return cell1.distance<cell2.distance;
}

void runForceDirectedPlacementAlgorithm(int noOfCells){
    int iterationLimit = 5;
    int iterationCount = 0;
    int abortCount = 0;
    int abortLimit = 0.35*noOfCells;
    bool endRipple = true;
    int targetX = 0;
    int targetY = 0;
    int connectedCellX = 0;
    int connectedCellY = 0;
    int connectedCellNumber = -1;
    cellWeightInfo currentCellWeightInfo;
    cellWeightInfo tempCellWeightInfo;
    cellInfo sourceCellInfo;
    cellInfo newSourceCellInfo;
    cellInfo targetCellInfo;
    cellInfo tempCellInfo;
    int counter = 0;
    int  cellCount = cellWeightVector.size();
    int cellConnectivitySize = 0;
    coOrdinates emptyCell;
    coOrdinates sourceCoOrdinates;
    //int counterTemp = 0;
    while(iterationCount<iterationLimit){
        currentCellWeightInfo = cellWeightVector.at(counter);
        endRipple = true;
        while(currentCellWeightInfo.blocked){
            counter++;
            if(counter<cellCount){
                currentCellWeightInfo = cellWeightVector.at(counter);
            }else{
                break;
            }
        }
        if(counter>=cellCount){
            for(int i=0;i<cellCount;i++){
                tempCellWeightInfo = cellWeightVector.at(i);
                tempCellWeightInfo.blocked = false;
                cellWeightVector.at(i) = tempCellWeightInfo;
                tempCellInfo = cellPositionVector.at(tempCellWeightInfo.cellNumber);
                tempCellInfo.blocked = false;
                cellPlacementMatrix[tempCellInfo.x][tempCellInfo.y] = tempCellInfo;
                cellPositionVector.at(tempCellInfo.cellNumber) = tempCellInfo;
            }
            counter = 0;
            currentCellWeightInfo = cellWeightVector.at(counter);
            iterationCount++;
        }
        if(!currentCellWeightInfo.blocked){
            sourceCellInfo = cellPositionVector.at(currentCellWeightInfo.cellNumber);
            endRipple = false;
        }
        while(endRipple == false){
            //printPlacementMatrix();
            cellConnectivitySize = cellConnectivityMatrix[currentCellWeightInfo.cellNumber].size();
            targetX = 0;
            targetY = 0;
            for(int i=0;i<cellConnectivitySize;i++){
                connectedCellNumber = cellConnectivityMatrix[currentCellWeightInfo.cellNumber][i];
                connectedCellX = cellPositionVector[connectedCellNumber].x;
                connectedCellY = cellPositionVector[connectedCellNumber].y;
                targetX = targetX + connectedCellX;
                targetY = targetY + connectedCellY;
            }
            if(cellConnectivitySize!=0){
                targetX = targetX/cellConnectivitySize;
                targetY = targetY/cellConnectivitySize;
                targetCellInfo = cellPlacementMatrix[targetX][targetY];
            }else{
                targetCellInfo = cellPlacementMatrix[targetX][targetY];
            }
            if(targetCellInfo.cellNumber == -1){//new position is vacant
                //placing the cell at target location
                targetCellInfo.cellNumber = currentCellWeightInfo.cellNumber;
                targetCellInfo.blocked = true;
                targetCellInfo.cellWeightVectorIndex = sourceCellInfo.cellWeightVectorIndex;
                currentCellWeightInfo.blocked = true;
                cellPlacementMatrix[targetCellInfo.x][targetCellInfo.y] = targetCellInfo;
                cellPositionVector.at(targetCellInfo.cellNumber) = targetCellInfo;
                cellWeightVector.at(targetCellInfo.cellWeightVectorIndex) = currentCellWeightInfo;
                //clearing source information
                if(!sourceCellInfo.blocked){
                    sourceCellInfo.blocked = false;
                    sourceCellInfo.cellNumber = -1;
                    sourceCellInfo.cellWeightVectorIndex = -1;
                    cellPlacementMatrix[sourceCellInfo.x][sourceCellInfo.y] = sourceCellInfo;
                }
                endRipple = true;
                abortCount = 0;
            }else if(!targetCellInfo.blocked && (targetCellInfo.x == sourceCellInfo.x && targetCellInfo.y == sourceCellInfo.y)){//same as current location
                sourceCellInfo.blocked = true;
                cellPlacementMatrix[sourceCellInfo.x][sourceCellInfo.y] = sourceCellInfo;
                cellPositionVector.at(sourceCellInfo.cellNumber) = sourceCellInfo;
                cellWeightVector.at(sourceCellInfo.cellWeightVectorIndex).blocked = true;
                endRipple = true;
                abortCount = 0;
            }else if(targetCellInfo.blocked){//target location is blocked
                sourceCoOrdinates.x = sourceCellInfo.x;
                sourceCoOrdinates.y = sourceCellInfo.y;
                emptyCell = findTargetEmptyLocation(sourceCoOrdinates);
                targetCellInfo = cellPlacementMatrix[emptyCell.x][emptyCell.y];
                targetCellInfo.cellNumber = currentCellWeightInfo.cellNumber;
                targetCellInfo.blocked = true;
                targetCellInfo.cellWeightVectorIndex = sourceCellInfo.cellWeightVectorIndex;
                currentCellWeightInfo.blocked = true;
                cellPlacementMatrix[targetCellInfo.x][targetCellInfo.y] = targetCellInfo;
                cellPositionVector.at(targetCellInfo.cellNumber) = targetCellInfo;
                cellWeightVector.at(targetCellInfo.cellWeightVectorIndex) = currentCellWeightInfo;
                //clearing source information
                if(!sourceCellInfo.blocked){
                    sourceCellInfo.blocked = false;
                    sourceCellInfo.cellNumber = -1;
                    sourceCellInfo.cellWeightVectorIndex = -1;
                    cellPlacementMatrix[sourceCellInfo.x][sourceCellInfo.y] = sourceCellInfo;
                }
                endRipple = true;
                abortCount++;
                if(abortCount>abortLimit){
                    for(int i=0;i<cellCount;i++){
                        tempCellWeightInfo = cellWeightVector.at(i);
                        tempCellWeightInfo.blocked = false;
                        cellWeightVector.at(i) = tempCellWeightInfo;
                        tempCellInfo = cellPositionVector.at(tempCellWeightInfo.cellNumber);
                        tempCellInfo.blocked = false;
                        cellPlacementMatrix[tempCellInfo.x][tempCellInfo.y] = tempCellInfo;
                        cellPositionVector.at(tempCellInfo.cellNumber) = tempCellInfo;
                    }
                    counter = 0;
                    iterationCount++;
                    break;
                }
            }else if(!targetCellInfo.blocked && targetCellInfo.cellNumber != -1){//target location is occupied
                newSourceCellInfo = targetCellInfo;
                //placing the cell at target location
                targetCellInfo.cellNumber = currentCellWeightInfo.cellNumber;
                targetCellInfo.blocked = true;
                newSourceCellInfo.blocked = true;
                targetCellInfo.cellWeightVectorIndex = sourceCellInfo.cellWeightVectorIndex;
                currentCellWeightInfo.blocked = true;
                cellPlacementMatrix[targetCellInfo.x][targetCellInfo.y] = targetCellInfo;
                cellPositionVector.at(targetCellInfo.cellNumber) = targetCellInfo;
                cellWeightVector.at(sourceCellInfo.cellWeightVectorIndex) = currentCellWeightInfo;
                //clearing source information
                if(!sourceCellInfo.blocked){
                    sourceCellInfo.blocked = false;
                    sourceCellInfo.cellNumber = -1;
                    sourceCellInfo.cellWeightVectorIndex = -1;
                    cellPlacementMatrix[sourceCellInfo.x][sourceCellInfo.y] = sourceCellInfo;
                }
                sourceCellInfo = newSourceCellInfo;
                currentCellWeightInfo = cellWeightVector.at(sourceCellInfo.cellWeightVectorIndex);
                endRipple = false;
                abortCount = 0;
            }
        }
    }
}

void prioritizeNets(){
    coOrdinates rect1;
    coOrdinates rect2;
    coOrdinates rect3;
    coOrdinates rect4;
    coOrdinates innerPoint1;
    coOrdinates innerPoint2;
    int netTerminalVectorSize = netTerminalVector.size();
    netInfo currentNetInfo;
    netInfo netToBeChecked;
    int priority = 0;
    bool isPointInsideRectangle = false;
    cellcoOrdinates sourceCellCoOrdinates;
    cellcoOrdinates targetCellCoOrdinates;
    cellcoOrdinates innersourceCellCoOrdinates;
    cellcoOrdinates innertargetCellCoOrdinates;
    for(int i=0;i<netTerminalVectorSize;i++){
        priority = 0;
        currentNetInfo = netTerminalVector.at(i);
        sourceCellCoOrdinates = cellCoOrdinateVector.at(currentNetInfo.sourceCellNumber-1);
        targetCellCoOrdinates = cellCoOrdinateVector.at(currentNetInfo.targetCellNumber-1);
        if(currentNetInfo.sourceCellTerminalNumber == 1){
            rect1 = sourceCellCoOrdinates.t1;
        }else if(currentNetInfo.sourceCellTerminalNumber == 2){
            rect1 = sourceCellCoOrdinates.t2;
        }else if(currentNetInfo.sourceCellTerminalNumber == 3){
            rect1 = sourceCellCoOrdinates.t3;
        }else if(currentNetInfo.sourceCellTerminalNumber == 4){
            rect1 = sourceCellCoOrdinates.t4;
        }
        if(currentNetInfo.targetCellTerminalNumber == 1){
            rect2 = targetCellCoOrdinates.t1;
        }else if(currentNetInfo.targetCellTerminalNumber == 2){
            rect2 = targetCellCoOrdinates.t2;
        }else if(currentNetInfo.targetCellTerminalNumber == 3){
            rect2 = targetCellCoOrdinates.t3;
        }else if(currentNetInfo.targetCellTerminalNumber == 4){
            rect2 = targetCellCoOrdinates.t4;
        }
        rect3.x = rect1.x;
        rect3.y = rect2.y;
        rect4.x = rect2.x;
        rect4.y = rect1.y;
        for(int j=0;j<netTerminalVectorSize;j++){
            if(i!=j){
                netToBeChecked = netTerminalVector.at(j);
                innersourceCellCoOrdinates = cellCoOrdinateVector.at(netToBeChecked.sourceCellNumber-1);
                innertargetCellCoOrdinates = cellCoOrdinateVector.at(netToBeChecked.targetCellNumber-1);
                if(netToBeChecked.sourceCellTerminalNumber == 1){
                    innerPoint1 = innersourceCellCoOrdinates.t1;
                }else if(netToBeChecked.sourceCellTerminalNumber == 2){
                    innerPoint1 = innersourceCellCoOrdinates.t2;
                }else if(netToBeChecked.sourceCellTerminalNumber == 3){
                    innerPoint1 = innersourceCellCoOrdinates.t3;
                }else if(netToBeChecked.sourceCellTerminalNumber == 4){
                    innerPoint1 = innersourceCellCoOrdinates.t4;
                }
                if(netToBeChecked.targetCellTerminalNumber == 1){
                    innerPoint2 = innertargetCellCoOrdinates.t1;
                }else if(netToBeChecked.targetCellTerminalNumber == 2){
                    innerPoint2 = innertargetCellCoOrdinates.t2;
                }else if(netToBeChecked.targetCellTerminalNumber == 3){
                    innerPoint2 = innertargetCellCoOrdinates.t3;
                }else if(netToBeChecked.targetCellTerminalNumber == 4){
                    innerPoint2 = innertargetCellCoOrdinates.t4;
                }
                isPointInsideRectangle = checkIfAPointIsInsideARectangle(rect1,rect2,rect3,rect4,innerPoint1);
                if(isPointInsideRectangle){
                    priority++;
                }else{
                    isPointInsideRectangle = checkIfAPointIsInsideARectangle(rect1,rect2,rect3,rect4,innerPoint2);
                    if(isPointInsideRectangle){
                        priority++;
                    }
                }
            }
        }
        currentNetInfo.netPriority = priority;
        netTerminalVector.at(i) = currentNetInfo;
    }
    sort(netTerminalVector.begin(),netTerminalVector.end(),compareNetInfo);
}
bool compareNetInfo(netInfo net1,netInfo net2){
    return net1.netPriority<net2.netPriority;
}
bool sortByNetNumber(netInfo net1,netInfo net2){
    return net1.netNumber<net2.netNumber;
}
bool checkIfAPointIsInsideARectangle(coOrdinates rect1,coOrdinates rect2,coOrdinates rect3,coOrdinates rect4,coOrdinates point){
    float trinagle1Area = areaOfTriangle(rect1,rect2,point);
    float trinagle2Area = areaOfTriangle(rect2,rect3,point);
    float trinagle3Area = areaOfTriangle(rect3,rect4,point);
    float trinagle4Area = areaOfTriangle(rect4,rect1,point);
    float rectArea = areaOfRect(rect1,rect2,rect3,rect4);
    float triangleAreaSum = (trinagle1Area+trinagle2Area+trinagle3Area+trinagle4Area);
    triangleAreaSum = ceil(triangleAreaSum);
    rectArea = ceil(rectArea);
    if(triangleAreaSum==rectArea){
        return true;
    }else{
        return false;
    }
}

float areaOfTriangle(coOrdinates a,coOrdinates b,coOrdinates p){
    float side1 = sqrt(pow(abs(a.y-b.y),2)+pow(abs(a.x-b.x),2));
    float side2 = sqrt(pow(abs(a.y-p.y),2)+pow(abs(a.x-p.x),2));
    float side3 = sqrt(pow(abs(b.y-p.y),2)+pow(abs(b.x-p.x),2));
    float semiPerimeter = (side1+side2+side3)/2;
    float area = sqrt(semiPerimeter*(semiPerimeter-side1)*(semiPerimeter-side2)*(semiPerimeter-side3));
    return area;
}

float areaOfRect(coOrdinates rect1,coOrdinates rect2,coOrdinates rect3,coOrdinates rect4){
    float side1 = sqrt(pow(abs(rect1.y-rect2.y),2)+pow(abs(rect1.x-rect2.x),2));
    float side2 = sqrt(pow(abs(rect2.y-rect3.y),2)+pow(abs(rect2.x-rect3.x),2));
    float area = side1*side2;
    return area;
}

coOrdinates getTerminalCoOrdinates(cellcoOrdinates cell,int terminalNumber){
    coOrdinates temp;
    if(terminalNumber == 1){
        return cell.t1;
    }else if(terminalNumber == 2){
        return cell.t2;
    }else if(terminalNumber == 3){
        return cell.t3;
    }else if(terminalNumber == 4){
        return cell.t4;
    }
    return temp;
}

int getFinalOrientation(int orientation){
    int finalOrientation = 0;
    if(orientation == 0){
        finalOrientation = 1;
    }else if(orientation == 1){
        finalOrientation = 2;
    }else if(orientation == 3){
        finalOrientation = 1;
    }
    return finalOrientation;
}
//0-top
//1-right
//2-down
//3-left
void flipACell(int cellNumber,int orientation){
    cellcoOrdinates currentCell = cellCoOrdinateVector.at(cellNumber);
    int currentOrientation = -1;
    coOrdinates referencePoint = currentCell.referencePoint;
    coOrdinates t1;
    coOrdinates t2;
    coOrdinates t3;
    coOrdinates t4;
    if(orientation == -1){
        currentOrientation = currentCell.orientation;
        if(currentOrientation == -1){
            currentCell.orientation = 0;
            currentOrientation = 0;
        }
    }else{
        currentOrientation = orientation;
    }
    if(currentOrientation == 0){
        currentCell.orientation = 0;
        t1.x = referencePoint.x;
        t1.y = referencePoint.y+1;
        currentCell.t1 = t1;
        t2.x = referencePoint.x;
        t2.y = referencePoint.y+4;
        currentCell.t2 = t2;
        t3.x = referencePoint.x+5;
        t3.y = referencePoint.y+1;
        currentCell.t3 = t3;
        t4.x = referencePoint.x+5;
        t4.y = referencePoint.y+4;
        currentCell.t4 = t4;
        /*chipGridMatrix[currentCell.t1.x-1][currentCell.t1.y].m1Blocked = 0;
        chipGridMatrix[currentCell.t1.x-1][currentCell.t1.y].m2Blocked = 0;
        chipGridMatrix[currentCell.t2.x-1][currentCell.t2.y].m1Blocked = 0;
        chipGridMatrix[currentCell.t2.x-1][currentCell.t2.y].m2Blocked = 0;
        chipGridMatrix[currentCell.t3.x+1][currentCell.t3.y].m1Blocked = 0;
        chipGridMatrix[currentCell.t3.x+1][currentCell.t3.y].m2Blocked = 0;
        chipGridMatrix[currentCell.t4.x+1][currentCell.t4.y].m1Blocked = 0;
        chipGridMatrix[currentCell.t4.x+1][currentCell.t4.y].m2Blocked = 0;*/
        for(int p=0;p<5;p++){
            chipGridMatrix[currentCell.t1.x-p][currentCell.t1.y].m1Blocked = 1;
            chipGridMatrix[currentCell.t1.x-p][currentCell.t1.y].m2Blocked = 1;
            chipGridMatrix[currentCell.t1.x-p][currentCell.t1.y].cellSurrounding = cellNumber;
            chipGridMatrix[currentCell.t1.x-p][currentCell.t1.y].cellTerminal = 1;
        }
        for(int p=0;p<3;p++){
            chipGridMatrix[currentCell.t2.x-p][currentCell.t2.y].m1Blocked = 1;
            chipGridMatrix[currentCell.t2.x-p][currentCell.t2.y].m2Blocked = 1;
            chipGridMatrix[currentCell.t2.x-p][currentCell.t2.y].cellSurrounding = cellNumber;
            chipGridMatrix[currentCell.t2.x-p][currentCell.t2.y].cellTerminal = 2;
        }
        for(int p=0;p<5;p++){
            chipGridMatrix[currentCell.t3.x+p][currentCell.t3.y].m1Blocked = 1;
            chipGridMatrix[currentCell.t3.x+p][currentCell.t3.y].m2Blocked = 1;
            chipGridMatrix[currentCell.t3.x+p][currentCell.t3.y].cellSurrounding = cellNumber;
            chipGridMatrix[currentCell.t3.x+p][currentCell.t3.y].cellTerminal = 3;
        }
        for(int p=0;p<3;p++){
            chipGridMatrix[currentCell.t4.x+p][currentCell.t4.y].m1Blocked = 1;
            chipGridMatrix[currentCell.t4.x+p][currentCell.t4.y].m2Blocked = 1;
            chipGridMatrix[currentCell.t4.x+p][currentCell.t4.y].cellSurrounding = cellNumber;
            chipGridMatrix[currentCell.t4.x+p][currentCell.t4.y].cellTerminal = 4;
        }
    }else if(currentOrientation == 1){
        currentCell.orientation = 1;
        t1.x = referencePoint.x+1;
        t1.y = referencePoint.y+5;
        currentCell.t1 = t1;
        t2.x = referencePoint.x+4;
        t2.y = referencePoint.y+5;
        currentCell.t2 = t2;
        t3.x = referencePoint.x+1;
        t3.y = referencePoint.y;
        currentCell.t3 = t3;
        t4.x = referencePoint.x+4;
        t4.y = referencePoint.y;
        currentCell.t4 = t4;
        /*chipGridMatrix[currentCell.t1.x][currentCell.t1.y+1].m1Blocked = 0;
        chipGridMatrix[currentCell.t1.x][currentCell.t1.y+1].m2Blocked = 0;
        chipGridMatrix[currentCell.t2.x][currentCell.t2.y+1].m1Blocked = 0;
        chipGridMatrix[currentCell.t2.x][currentCell.t2.y+1].m2Blocked = 0;
        chipGridMatrix[currentCell.t3.x][currentCell.t3.y-1].m1Blocked = 0;
        chipGridMatrix[currentCell.t3.x][currentCell.t3.y-1].m2Blocked = 0;
        chipGridMatrix[currentCell.t4.x][currentCell.t4.y-1].m1Blocked = 0;
        chipGridMatrix[currentCell.t4.x][currentCell.t4.y-1].m2Blocked = 0;*/
        for(int p=0;p<5;p++){
            chipGridMatrix[currentCell.t1.x][currentCell.t1.y+p].m1Blocked = 1;
            chipGridMatrix[currentCell.t1.x][currentCell.t1.y+p].m2Blocked = 1;
            chipGridMatrix[currentCell.t1.x][currentCell.t1.y+p].cellSurrounding = cellNumber;
            chipGridMatrix[currentCell.t1.x][currentCell.t1.y+p].cellTerminal = 1;
        }
        for(int p=0;p<3;p++){
            chipGridMatrix[currentCell.t2.x][currentCell.t2.y+p].m1Blocked = 1;
            chipGridMatrix[currentCell.t2.x][currentCell.t2.y+p].m2Blocked = 1;
            chipGridMatrix[currentCell.t2.x][currentCell.t2.y+p].cellSurrounding = cellNumber;
            chipGridMatrix[currentCell.t2.x][currentCell.t2.y+p].cellTerminal = 2;
        }
        for(int p=0;p<5;p++){
            chipGridMatrix[currentCell.t3.x][currentCell.t3.y-p].m1Blocked = 1;
            chipGridMatrix[currentCell.t3.x][currentCell.t3.y-p].m2Blocked = 1;
            chipGridMatrix[currentCell.t3.x][currentCell.t3.y-p].cellSurrounding = cellNumber;
            chipGridMatrix[currentCell.t3.x][currentCell.t3.y-p].cellTerminal = 3;
        }
        for(int p=0;p<3;p++){
            chipGridMatrix[currentCell.t4.x][currentCell.t4.y-p].m1Blocked = 1;
            chipGridMatrix[currentCell.t4.x][currentCell.t4.y-p].m2Blocked = 1;
            chipGridMatrix[currentCell.t4.x][currentCell.t4.y-p].cellSurrounding = cellNumber;
            chipGridMatrix[currentCell.t4.x][currentCell.t4.y-p].cellTerminal = 4;
        }
    }else if(currentOrientation == 2){
        currentCell.orientation = 2;
        t1.x = referencePoint.x+5;
        t1.y = referencePoint.y+4;
        currentCell.t1 = t1;
        t2.x = referencePoint.x+5;
        t2.y = referencePoint.y+1;
        currentCell.t2 = t2;
        t3.x = referencePoint.x;
        t3.y = referencePoint.y+4;
        currentCell.t3 = t3;
        t4.x = referencePoint.x;
        t4.y = referencePoint.y+1;
        currentCell.t4 = t4;
        /*chipGridMatrix[currentCell.t1.x+1][currentCell.t1.y].m1Blocked = 0;
        chipGridMatrix[currentCell.t1.x+1][currentCell.t1.y].m2Blocked = 0;
        chipGridMatrix[currentCell.t2.x+1][currentCell.t2.y].m1Blocked = 0;
        chipGridMatrix[currentCell.t2.x+1][currentCell.t2.y].m2Blocked = 0;
        chipGridMatrix[currentCell.t3.x-1][currentCell.t3.y].m1Blocked = 0;
        chipGridMatrix[currentCell.t3.x-1][currentCell.t3.y].m2Blocked = 0;
        chipGridMatrix[currentCell.t4.x-1][currentCell.t4.y].m1Blocked = 0;
        chipGridMatrix[currentCell.t4.x-1][currentCell.t4.y].m2Blocked = 0;*/
        for(int p=0;p<5;p++){
            chipGridMatrix[currentCell.t1.x+p][currentCell.t1.y].m1Blocked = 1;
            chipGridMatrix[currentCell.t1.x+p][currentCell.t1.y].m2Blocked = 1;
            chipGridMatrix[currentCell.t1.x+p][currentCell.t1.y].cellSurrounding = cellNumber;
            chipGridMatrix[currentCell.t1.x+p][currentCell.t1.y].cellTerminal = 1;
        }
        for(int p=0;p<3;p++){
            chipGridMatrix[currentCell.t2.x+p][currentCell.t2.y].m1Blocked = 1;
            chipGridMatrix[currentCell.t2.x+p][currentCell.t2.y].m2Blocked = 1;
            chipGridMatrix[currentCell.t2.x+p][currentCell.t2.y].cellSurrounding = cellNumber;
            chipGridMatrix[currentCell.t2.x+p][currentCell.t2.y].cellTerminal = 2;
        }
        for(int p=0;p<5;p++){
            chipGridMatrix[currentCell.t3.x-p][currentCell.t3.y].m1Blocked = 1;
            chipGridMatrix[currentCell.t3.x-p][currentCell.t3.y].m2Blocked = 1;
            chipGridMatrix[currentCell.t3.x-p][currentCell.t3.y].cellSurrounding = cellNumber;
            chipGridMatrix[currentCell.t3.x-p][currentCell.t3.y].cellTerminal = 3;
        }
        for(int p=0;p<3;p++){
            chipGridMatrix[currentCell.t4.x-p][currentCell.t4.y].m1Blocked = 1;
            chipGridMatrix[currentCell.t4.x-p][currentCell.t4.y].m2Blocked = 1;
            chipGridMatrix[currentCell.t4.x-p][currentCell.t4.y].cellSurrounding = cellNumber;
            chipGridMatrix[currentCell.t4.x-p][currentCell.t4.y].cellTerminal = 4;
        }

    }else if(currentOrientation == 3){
        currentCell.orientation = 3;
        t1.x = referencePoint.x+4;
        t1.y = referencePoint.y;
        currentCell.t1 = t1;
        t2.x = referencePoint.x+1;
        t2.y = referencePoint.y;
        currentCell.t2 = t2;
        t3.x = referencePoint.x+5;
        t3.y = referencePoint.y+4;
        currentCell.t3 = t3;
        t4.x = referencePoint.x+5;
        t4.y = referencePoint.y+4;
        currentCell.t4 = t4;
        /*chipGridMatrix[currentCell.t1.x][currentCell.t1.y-1].m1Blocked = 0;
        chipGridMatrix[currentCell.t1.x][currentCell.t1.y-1].m2Blocked = 0;
        chipGridMatrix[currentCell.t2.x][currentCell.t2.y-1].m1Blocked = 0;
        chipGridMatrix[currentCell.t2.x][currentCell.t2.y-1].m2Blocked = 0;
        chipGridMatrix[currentCell.t3.x][currentCell.t3.y-1].m1Blocked = 0;
        chipGridMatrix[currentCell.t3.x][currentCell.t3.y+1].m2Blocked = 0;
        chipGridMatrix[currentCell.t4.x][currentCell.t4.y+1].m1Blocked = 0;
        chipGridMatrix[currentCell.t4.x][currentCell.t4.y+1].m2Blocked = 0;*/
        for(int p=0;p<5;p++){
            chipGridMatrix[currentCell.t1.x][currentCell.t1.y-p].m1Blocked = 1;
            chipGridMatrix[currentCell.t1.x][currentCell.t1.y-p].m2Blocked = 1;
            chipGridMatrix[currentCell.t1.x][currentCell.t1.y-p].cellSurrounding = cellNumber;
            chipGridMatrix[currentCell.t1.x][currentCell.t1.y-p].cellTerminal = 1;
        }
        for(int p=0;p<3;p++){
            chipGridMatrix[currentCell.t2.x][currentCell.t2.y-p].m1Blocked = 1;
            chipGridMatrix[currentCell.t2.x][currentCell.t2.y-p].m2Blocked = 1;
            chipGridMatrix[currentCell.t2.x][currentCell.t2.y-p].cellSurrounding = cellNumber;
            chipGridMatrix[currentCell.t2.x][currentCell.t2.y-p].cellTerminal = 2;
        }
        for(int p=0;p<5;p++){
            chipGridMatrix[currentCell.t3.x][currentCell.t3.y+p].m1Blocked = 1;
            chipGridMatrix[currentCell.t3.x][currentCell.t3.y+p].m2Blocked = 1;
            chipGridMatrix[currentCell.t3.x][currentCell.t3.y+p].cellSurrounding = cellNumber;
            chipGridMatrix[currentCell.t3.x][currentCell.t3.y+p].cellTerminal = 3;
        }
        for(int p=0;p<3;p++){
            chipGridMatrix[currentCell.t4.x][currentCell.t4.y+p].m1Blocked = 1;
            chipGridMatrix[currentCell.t4.x][currentCell.t4.y+p].m2Blocked = 1;
            chipGridMatrix[currentCell.t4.x][currentCell.t4.y+p].cellSurrounding = cellNumber;
            chipGridMatrix[currentCell.t4.x][currentCell.t4.y+p].cellTerminal = 4;
        }
    }
    cellCoOrdinateVector.at(cellNumber) = currentCell;
}

void convertPlacementMatrixToChipGridMatrix(int noOfCells,int cellToBeFlipped,int rowSpacing,int columnSpacing){
    int placementRows = cellPlacementMatrix.size();
    int placementColumns = cellPlacementMatrix.size();
    int chipGridMatrixRows = placementRows*(11+rowSpacing)+rowSpacing;
    int chipGridMatrixColumns = placementColumns*(11+columnSpacing)+columnSpacing;
    chipGridMatrix = vector< vector<chipUnitBlock> >(chipGridMatrixRows,vector<chipUnitBlock>(chipGridMatrixColumns));
    chipUnitBlock currentChipUnitBlock;
    cellcoOrdinates flippedCellCoOrdinates;
    int newOrientation = -1;
    coOrdinates point;
    cellCoOrdinateVector.resize(noOfCells);
    int cellNumber = -1;
    for(int i=0;i<chipGridMatrixRows;i++){
        for(int j=0;j<chipGridMatrixColumns;j++){
            currentChipUnitBlock.m1Blocked = 0;
            currentChipUnitBlock.m2Blocked = 0;
            currentChipUnitBlock.cellNumber = -1;
            currentChipUnitBlock.cellSurrounding = -1;
            currentChipUnitBlock.m1NetSurrounding = 0;
            currentChipUnitBlock.m2NetSurrounding = 0;
            currentChipUnitBlock.cellTerminal = -1;
            currentChipUnitBlock.m1NetNumber= 0;
            currentChipUnitBlock.m2NetNumber = 0;
            currentChipUnitBlock.m1LeeNumber = -1;
            currentChipUnitBlock.m2LeeNumber = -1;
            currentChipUnitBlock.via = 0;
            currentChipUnitBlock.x = i;
            currentChipUnitBlock.y = j;
            chipGridMatrix[i][j] = currentChipUnitBlock;
        }
    }
    for(int i=0;i<placementRows;i++){
        for(int j=0;j<placementColumns;j++){
            point.x = i*6+rowSpacing*(i+1);
            point.y = j*6+columnSpacing*(j+1);
            cellNumber = cellPlacementMatrix[i][j].cellNumber;
            if(cellNumber!=-1){
                cellCoOrdinateVector.at(cellNumber).referencePoint = point;
                for(int u=point.x-1;u<point.x+7;u++){
                    for(int v=point.y-1;v<point.y+7;v++){
                        if(!(u==point.x-1 || u==point.x+6 || v==point.y-1 ||v==point.y+6)){
                            chipGridMatrix[u][v].cellNumber = cellNumber;
                        }
                        chipGridMatrix[u][v].m1Blocked = 1;
                        chipGridMatrix[u][v].m2Blocked = 1;
                    }
                }
                flipACell(cellNumber,-1);
                if(cellNumber == cellToBeFlipped){
                    flippedCellCoOrdinates = cellCoOrdinateVector.at(cellToBeFlipped);
                    newOrientation = getFinalOrientation(flippedCellCoOrdinates.orientation);
                    flipACell(cellNumber,newOrientation);
                }
            }
        }
    }
}

bool compareCellInfo(cellWeightInfo node1,cellWeightInfo node2){
    return node1.cellWeight>node2.cellWeight;
}


bool isCompletelyBlocked(chipUnitBlock tempChipUnitBlock){
    return (tempChipUnitBlock.m1Blocked == 1 && tempChipUnitBlock.m2Blocked == 1);
}

bool isBlockInSurrounding(chipUnitBlock currentBlock,chipUnitBlock source,chipUnitBlock target){
    return (currentBlock.cellSurrounding == source.cellNumber && currentBlock.cellTerminal == source.cellTerminal) || (currentBlock.cellSurrounding == target.cellNumber && currentBlock.cellTerminal == target.cellTerminal);
}
void clearNetSurroundingNumber(int metalType,int netNumber,chipUnitBlock currentBlock){
    if(metalType == 1 && currentBlock.m1NetSurrounding == netNumber){
        currentBlock.m1Blocked = 0;
        currentBlock.m1NetSurrounding = 0;
    }else if(metalType == 2 && currentBlock.m2NetSurrounding == netNumber){
        currentBlock.m2Blocked = 0;
        currentBlock.m2NetSurrounding = 0;
    }
    chipGridMatrix[currentBlock.x][currentBlock.y] = currentBlock;
}

void updateNetSurroundingNumber(int metalType,int netNumber,chipUnitBlock currentBlock){
    if(metalType == 1 && currentBlock.m1NetSurrounding == 0){
        currentBlock.m1Blocked = 1;
        currentBlock.m1NetSurrounding = netNumber;
    }else if(metalType == 2 && currentBlock.m2NetSurrounding == 0){
        currentBlock.m2Blocked = 1;
        currentBlock.m2NetSurrounding = netNumber;
    }
    chipGridMatrix[currentBlock.x][currentBlock.y] = currentBlock;
}

int performLeeAlgorithmForTheNets(int noOfNets,int &failedMethod, int metalType){
    netInfo currentNetInfo;
    int sourceCellNumber = -1;
    int targetCellNumber = -1;
    chipUnitBlock source;
    chipUnitBlock target;
    int sourceX = -1;
    int sourceY = -1;
    int targetX = -1;
    int targetY = -1;
    int netNumber = 0;
    cellcoOrdinates sourceCellCoOrdinates;
    coOrdinates sourceTerminalCoOrdinates;
    coOrdinates targetTerminalCoOrdinates;
    int counter = 0;
    for(int i=0;i<noOfNets;i++){
        currentNetInfo = netTerminalVector.at(i);
        if(currentNetInfo.routingMetal == -1){
            netNumber = currentNetInfo.netNumber;
            //if(netNumber == 1 || netNumber == 3){
            sourceCellNumber = currentNetInfo.sourceCellNumber;
            targetCellNumber = currentNetInfo.targetCellNumber;
            sourceCellCoOrdinates = cellCoOrdinateVector.at(sourceCellNumber-1);
            sourceTerminalCoOrdinates = getTerminalCoOrdinates(cellCoOrdinateVector.at(sourceCellNumber-1),currentNetInfo.sourceCellTerminalNumber);
            source = chipGridMatrix[sourceTerminalCoOrdinates.x][sourceTerminalCoOrdinates.y];
            sourceX = source.x;
            sourceY = source.y;
            targetTerminalCoOrdinates = getTerminalCoOrdinates(cellCoOrdinateVector.at(targetCellNumber-1),currentNetInfo.targetCellTerminalNumber);
            target = chipGridMatrix[targetTerminalCoOrdinates.x][targetTerminalCoOrdinates.y];
            targetX = target.x;
            targetY = target.y;
            failedMethod = generateLeeNumbersFromSourceToTerminal(source,target,metalType);//-1 failed in this method
            if(failedMethod == 0){
                failedMethod = backTraceNetPath(netNumber,chipGridMatrix[sourceX][sourceY],chipGridMatrix[targetX][targetY]);//-2 failed in this method
            }
            clearLeeNumbersFromSourceToTerminal(netNumber,chipGridMatrix[sourceX][sourceY],chipGridMatrix[targetX][targetY]);
            if(failedMethod == -2 || failedMethod == -1){
            }else{
                netTerminalVector.at(i).routingMetal = metalType;
                counter++;
                //cout<<"nets routed till now"<<i+1<<endl;
            }
        }
    //}//
    }
    cout<<"# of nets routed:"<< counter<<endl;
    return -1;
}
void printChipGridMatrix(){
    cout<<"ChipGrid Matrix:"<<endl;
    int rowSize = chipGridMatrix.size();
    int columnSize = -1;
    for(int i=59;i<120;i++){
        columnSize = chipGridMatrix[i].size();
        for(int j=45;j<80;j++){
            if(i==89 && j==70){
                cout<<" S"<<" ";
            }else if(i==64 && j==49){
                cout<<" T"<<" ";
            }else{
                //if(chipGridMatrix[i][j].m1LeeNumber== -1 && chipGridMatrix[i][j].m2LeeNumber== -1){
                //    cout<<". ";
                //}else
                if(chipGridMatrix[i][j].m1LeeNumber== -1){
                    cout<<chipGridMatrix[i][j].m2LeeNumber<<" ";
                }else if(chipGridMatrix[i][j].m2LeeNumber== -1){
                   cout<<chipGridMatrix[i][j].m1LeeNumber<<" ";
                }
            }
            //cout<<"("<<chipGridMatrix[i][j].m1Blocked<<","<<chipGridMatrix[i][j].m2Blocked<<","<<chipGridMatrix[i][j].cellNumber<<","<<chipGridMatrix[i][j].netNumber<<","<<chipGridMatrix[i][j].leeNumber<<")   ";
            //cout<<"("<<chipGridMatrix[i][j].cellNumber<<","<<chipGridMatrix[i][j].m1LeeNumber<<") ";
        }
        cout<<endl;
    }
}

void printChipGridMatrixElement(int i,int j){
    chipUnitBlock c = chipGridMatrix[i][j];
    cout<<"cellNumber "<<c.cellNumber<<endl;
    cout<<"cellSurrounding "<<c.cellSurrounding<<endl;
    cout<<"cellTerminal "<<c.cellTerminal<<endl;
    cout<<"m1Blocked "<<c.m1Blocked<<endl;
    cout<<"m2Blocked "<<c.m2Blocked<<endl;
    cout<<"m1LeeNumber "<<c.m1LeeNumber<<endl;
    cout<<"m1NetNumber "<<c.m1NetNumber<<endl;
    cout<<"m1NetSurrounding "<<c.m1NetSurrounding<<endl;
    cout<<"m2LeeNumber "<<c.m2LeeNumber<<endl;
    cout<<"m2NetNumber "<<c.m2NetNumber<<endl;
    cout<<"m2NetSurrounding "<<c.m2NetSurrounding<<endl;
    cout<<"m2NetSurrounding "<<c.via<<endl;
    cout<<"x "<<c.x<<endl;
    cout<<"y "<<c.y<<endl;
}

void printPlacementMatrix(){
    cout<<"Placement Matrix:"<<endl;
    int rowSize = cellPlacementMatrix.size();
    int columnSize = cellPlacementMatrix[0].size();
    cellInfo sourceCellInfo;
    for(int i=0;i<rowSize;i++){
        for(int j=0;j<columnSize;j++){
            sourceCellInfo = cellPlacementMatrix[i][j];
            //cout<<sourceCellInfo.blocked<<","<<sourceCellInfo.cellNumber<<","<<sourceCellInfo.cellWeightVectorIndex<<"     ";
            cout<<sourceCellInfo.cellNumber<<"   ";
        }
        cout<<endl;
    }
}

int generateLeeNumbersFromSourceToTerminal(chipUnitBlock source,chipUnitBlock target,int metalType){
    int sourceX = source.x;
    int sourceY = source.y;
    int sourceLeeNumber = 0;
    chipGridMatrix[sourceX][sourceY].m1LeeNumber = sourceLeeNumber;
    chipGridMatrix[sourceX][sourceY].m2LeeNumber = sourceLeeNumber;
    source.m1LeeNumber = sourceLeeNumber;
    source.m2LeeNumber = sourceLeeNumber;
    int currentX = -1;
    int currentY = -1;
    chipUnitBlock newSource;
    unVisitedQueue.push(source);
    int targetReached = -1;
    while(!unVisitedQueue.empty()){
        newSource = unVisitedQueue.front();
        sourceX = newSource.x;
        sourceY = newSource.y;
        if(newSource.m1LeeNumber == -1){
            sourceLeeNumber = newSource.m2LeeNumber;
        }else{
            sourceLeeNumber = newSource.m1LeeNumber;
        }
        if(newSource.x == target.x && newSource.y == target.y){
            targetReached = 0;
            break;
        }
        currentX = sourceX+1;
        currentY = sourceY;
        if(isPointInGridBounds(currentX,currentY) && chipGridMatrix[currentX][currentY].via==0 &&(!isCompletelyBlocked(chipGridMatrix[currentX][currentY])|| (currentX == target.x && currentY == target.y) || (isBlockInSurrounding(chipGridMatrix[currentX][currentY],source,target)))
           && chipGridMatrix[currentX][currentY].m1LeeNumber == -1 && chipGridMatrix[currentX][currentY].m2LeeNumber == -1){
            /*if(isBlockInSurrounding(chipGridMatrix[currentX][currentY],source,target)){
                chipGridMatrix[currentX][currentY].m1LeeNumber = sourceLeeNumber+1;
                chipGridMatrix[currentX][currentY].m2LeeNumber = sourceLeeNumber+1;
            }
            if(chipGridMatrix[currentX][currentY].m1Blocked == 0){
                chipGridMatrix[currentX][currentY].m1LeeNumber = sourceLeeNumber+1;
            }
            if(chipGridMatrix[currentX][currentY].m2Blocked == 0){
                chipGridMatrix[currentX][currentY].m2LeeNumber = sourceLeeNumber+1;
            }*/
            if(chipGridMatrix[currentX][currentY].m1Blocked == 0 ||(chipGridMatrix[currentX][currentY].m1NetSurrounding == 0 && isBlockInSurrounding(chipGridMatrix[currentX][currentY],source,target))){
                chipGridMatrix[currentX][currentY].m1LeeNumber = sourceLeeNumber+1;
            }else if(chipGridMatrix[currentX][currentY].m2Blocked == 0 ||(chipGridMatrix[currentX][currentY].m2NetSurrounding == 0 && isBlockInSurrounding(chipGridMatrix[currentX][currentY],source,target))){
                chipGridMatrix[currentX][currentY].m2LeeNumber = sourceLeeNumber+1;
            }
            //chipGridMatrix[currentX][currentY].leeNumber = sourceLeeNumber+1;
            unVisitedQueue.push(chipGridMatrix[currentX][currentY]);
        }
        currentX = sourceX-1;
        currentY = sourceY;
        if(isPointInGridBounds(currentX,currentY) && chipGridMatrix[currentX][currentY].via==0 && (!isCompletelyBlocked(chipGridMatrix[currentX][currentY]) || (currentX == target.x && currentY == target.y) || (isBlockInSurrounding(chipGridMatrix[currentX][currentY],source,target)))
           && chipGridMatrix[currentX][currentY].m1LeeNumber == -1 && chipGridMatrix[currentX][currentY].m2LeeNumber == -1){

            /*if(isBlockInSurrounding(chipGridMatrix[currentX][currentY],source,target)){
                chipGridMatrix[currentX][currentY].m1LeeNumber = sourceLeeNumber+1;
                chipGridMatrix[currentX][currentY].m2LeeNumber = sourceLeeNumber+1;
            }
            if(chipGridMatrix[currentX][currentY].m1Blocked == 0){
                chipGridMatrix[currentX][currentY].m1LeeNumber = sourceLeeNumber+1;
            }
            if(chipGridMatrix[currentX][currentY].m2Blocked == 0){
                chipGridMatrix[currentX][currentY].m2LeeNumber = sourceLeeNumber+1;
            }*/
            if(chipGridMatrix[currentX][currentY].m1Blocked == 0 || (chipGridMatrix[currentX][currentY].m1NetSurrounding == 0 && isBlockInSurrounding(chipGridMatrix[currentX][currentY],source,target))){
                chipGridMatrix[currentX][currentY].m1LeeNumber = sourceLeeNumber+1;
            }else if(chipGridMatrix[currentX][currentY].m2Blocked == 0 ||(chipGridMatrix[currentX][currentY].m2NetSurrounding == 0 && isBlockInSurrounding(chipGridMatrix[currentX][currentY],source,target))){
                chipGridMatrix[currentX][currentY].m2LeeNumber = sourceLeeNumber+1;
            }
            //chipGridMatrix[currentX][currentY].leeNumber = sourceLeeNumber+1;
            unVisitedQueue.push(chipGridMatrix[currentX][currentY]);
        }
        currentX = sourceX;
        currentY = sourceY+1;
        if(isPointInGridBounds(currentX,currentY) && chipGridMatrix[currentX][currentY].via==0 && (!isCompletelyBlocked(chipGridMatrix[currentX][currentY]) || (currentX == target.x && currentY == target.y) || (isBlockInSurrounding(chipGridMatrix[currentX][currentY],source,target)))
           && chipGridMatrix[currentX][currentY].m1LeeNumber == -1 && chipGridMatrix[currentX][currentY].m2LeeNumber == -1){

            /*if(isBlockInSurrounding(chipGridMatrix[currentX][currentY],source,target)){
                chipGridMatrix[currentX][currentY].m1LeeNumber = sourceLeeNumber+1;
                chipGridMatrix[currentX][currentY].m2LeeNumber = sourceLeeNumber+1;
            }
            if(chipGridMatrix[currentX][currentY].m1Blocked == 0){
                chipGridMatrix[currentX][currentY].m1LeeNumber = sourceLeeNumber+1;
            }
            if(chipGridMatrix[currentX][currentY].m2Blocked == 0){
                chipGridMatrix[currentX][currentY].m2LeeNumber = sourceLeeNumber+1;
            }*/
            if(chipGridMatrix[currentX][currentY].m1Blocked == 0 || (chipGridMatrix[currentX][currentY].m1NetSurrounding == 0 && isBlockInSurrounding(chipGridMatrix[currentX][currentY],source,target))){
                chipGridMatrix[currentX][currentY].m1LeeNumber = sourceLeeNumber+1;
            }else if(chipGridMatrix[currentX][currentY].m2Blocked == 0 || (chipGridMatrix[currentX][currentY].m2NetSurrounding == 0 && isBlockInSurrounding(chipGridMatrix[currentX][currentY],source,target))){
                chipGridMatrix[currentX][currentY].m2LeeNumber = sourceLeeNumber+1;
            }
            //chipGridMatrix[currentX][currentY].leeNumber = sourceLeeNumber+1;
            unVisitedQueue.push(chipGridMatrix[currentX][currentY]);
        }
        currentX = sourceX;
        currentY = sourceY-1;
        if(isPointInGridBounds(currentX,currentY) && chipGridMatrix[currentX][currentY].via==0 && (!isCompletelyBlocked(chipGridMatrix[currentX][currentY]) || (currentX == target.x && currentY == target.y) || (isBlockInSurrounding(chipGridMatrix[currentX][currentY],source,target)))
           && chipGridMatrix[currentX][currentY].m1LeeNumber == -1 && chipGridMatrix[currentX][currentY].m2LeeNumber == -1){

            /*if(isBlockInSurrounding(chipGridMatrix[currentX][currentY],source,target)){
                chipGridMatrix[currentX][currentY].m1LeeNumber = sourceLeeNumber+1;
                chipGridMatrix[currentX][currentY].m2LeeNumber = sourceLeeNumber+1;
            }
            if(chipGridMatrix[currentX][currentY].m1Blocked == 0){
                chipGridMatrix[currentX][currentY].m1LeeNumber = sourceLeeNumber+1;
            }
            if(chipGridMatrix[currentX][currentY].m2Blocked == 0){
                chipGridMatrix[currentX][currentY].m2LeeNumber = sourceLeeNumber+1;
            }*/
            if(chipGridMatrix[currentX][currentY].m1Blocked == 0 || (chipGridMatrix[currentX][currentY].m1NetSurrounding == 0 && isBlockInSurrounding(chipGridMatrix[currentX][currentY],source,target))){
                chipGridMatrix[currentX][currentY].m1LeeNumber = sourceLeeNumber+1;
            }else if(chipGridMatrix[currentX][currentY].m2Blocked == 0 || (chipGridMatrix[currentX][currentY].m2NetSurrounding == 0 && isBlockInSurrounding(chipGridMatrix[currentX][currentY],source,target))){
                chipGridMatrix[currentX][currentY].m2LeeNumber = sourceLeeNumber+1;
            }
            //chipGridMatrix[currentX][currentY].leeNumber = sourceLeeNumber+1;
            unVisitedQueue.push(chipGridMatrix[currentX][currentY]);
        }
        unVisitedQueue.pop();
    }
    while(!unVisitedQueue.empty()){
        unVisitedQueue.pop();
    }
    return targetReached;
}
bool isPointInGridBounds(int x, int y){
    int chipGridMatrixSize = chipGridMatrix.size();
    return !(x>=chipGridMatrixSize || x<0 || y>= chipGridMatrixSize || y<0);
}

void blockCellsAroundTheBendForMetal1(int netNumber,chipUnitBlock bendBlock){
    updateNetSurroundingNumber(1,netNumber,chipGridMatrix[bendBlock.x][bendBlock.y]);
    if(isPointInGridBounds(bendBlock.x+1,bendBlock.y)){
        updateNetSurroundingNumber(1,netNumber,chipGridMatrix[bendBlock.x+1][bendBlock.y]);
    }
    if(isPointInGridBounds(bendBlock.x-1,bendBlock.y)){
        updateNetSurroundingNumber(1,netNumber,chipGridMatrix[bendBlock.x-1][bendBlock.y]);

    }
    if(isPointInGridBounds(bendBlock.x,bendBlock.y-1)){
        updateNetSurroundingNumber(1,netNumber,chipGridMatrix[bendBlock.x][bendBlock.y-1]);
    }
    if(isPointInGridBounds(bendBlock.x,bendBlock.y+1)){
        updateNetSurroundingNumber(1,netNumber,chipGridMatrix[bendBlock.x][bendBlock.y+1]);
    }
    if(isPointInGridBounds(bendBlock.x+1,bendBlock.y-1)){
        updateNetSurroundingNumber(1,netNumber,chipGridMatrix[bendBlock.x+1][bendBlock.y-1]);
    }
    if(isPointInGridBounds(bendBlock.x+1,bendBlock.y+1)){
        updateNetSurroundingNumber(1,netNumber,chipGridMatrix[bendBlock.x+1][bendBlock.y+1]);
    }
    if(isPointInGridBounds(bendBlock.x-1,bendBlock.y+1)){
        updateNetSurroundingNumber(1,netNumber,chipGridMatrix[bendBlock.x-1][bendBlock.y+1]);
    }
    if(isPointInGridBounds(bendBlock.x-1,bendBlock.y-1)){
        chipGridMatrix[bendBlock.x-1][bendBlock.y-1].m1Blocked = 1;
    }
}

void unBlockCellsAroundTheBendForMetal1(int netNumber,chipUnitBlock bendBlock,chipUnitBlock source,chipUnitBlock target){
    clearNetSurroundingNumber(1,netNumber,chipGridMatrix[bendBlock.x][bendBlock.y]);
    if(isPointInGridBounds(bendBlock.x+1,bendBlock.y) && !isBlockInSurrounding(chipGridMatrix[bendBlock.x+1][bendBlock.y],source,target)){
        clearNetSurroundingNumber(1,netNumber,chipGridMatrix[bendBlock.x+1][bendBlock.y]);
    }
    if(isPointInGridBounds(bendBlock.x-1,bendBlock.y) && !isBlockInSurrounding(chipGridMatrix[bendBlock.x-1][bendBlock.y],source,target)){
        clearNetSurroundingNumber(1,netNumber,chipGridMatrix[bendBlock.x-1][bendBlock.y]);

    }
    if(isPointInGridBounds(bendBlock.x,bendBlock.y-1) && !isBlockInSurrounding(chipGridMatrix[bendBlock.x][bendBlock.y-1],source,target)){
        clearNetSurroundingNumber(1,netNumber,chipGridMatrix[bendBlock.x][bendBlock.y-1]);
    }
    if(isPointInGridBounds(bendBlock.x,bendBlock.y+1) && !isBlockInSurrounding(chipGridMatrix[bendBlock.x][bendBlock.y+1],source,target)){
        clearNetSurroundingNumber(1,netNumber,chipGridMatrix[bendBlock.x][bendBlock.y+1]);
    }
    if(isPointInGridBounds(bendBlock.x+1,bendBlock.y-1) && !isBlockInSurrounding(chipGridMatrix[bendBlock.x+1][bendBlock.y-1],source,target)){
        clearNetSurroundingNumber(1,netNumber,chipGridMatrix[bendBlock.x+1][bendBlock.y-1]);
    }
    if(isPointInGridBounds(bendBlock.x+1,bendBlock.y+1) && !isBlockInSurrounding(chipGridMatrix[bendBlock.x+1][bendBlock.y+1],source,target)){
        clearNetSurroundingNumber(1,netNumber,chipGridMatrix[bendBlock.x+1][bendBlock.y+1]);
    }
    if(isPointInGridBounds(bendBlock.x-1,bendBlock.y+1) && !isBlockInSurrounding(chipGridMatrix[bendBlock.x-1][bendBlock.y+1],source,target)){
        clearNetSurroundingNumber(1,netNumber,chipGridMatrix[bendBlock.x-1][bendBlock.y+1]);
    }
    if(isPointInGridBounds(bendBlock.x-1,bendBlock.y-1) && !isBlockInSurrounding(chipGridMatrix[bendBlock.x-1][bendBlock.y-1],source,target)){
        clearNetSurroundingNumber(1,netNumber,chipGridMatrix[bendBlock.x-1][bendBlock.y-1]);
    }
}

void blockCellsAroundTheBendForMetal2(int netNumber,chipUnitBlock bendBlock){
    updateNetSurroundingNumber(2,netNumber,chipGridMatrix[bendBlock.x][bendBlock.y]);
    if(isPointInGridBounds(bendBlock.x+1,bendBlock.y)){
        updateNetSurroundingNumber(2,netNumber,chipGridMatrix[bendBlock.x+1][bendBlock.y]);
    }
    if(isPointInGridBounds(bendBlock.x-1,bendBlock.y)){
        updateNetSurroundingNumber(2,netNumber,chipGridMatrix[bendBlock.x-1][bendBlock.y]);

    }
    if(isPointInGridBounds(bendBlock.x,bendBlock.y-1)){
        updateNetSurroundingNumber(2,netNumber,chipGridMatrix[bendBlock.x][bendBlock.y-1]);
    }
    if(isPointInGridBounds(bendBlock.x,bendBlock.y+1)){
        updateNetSurroundingNumber(2,netNumber,chipGridMatrix[bendBlock.x][bendBlock.y+1]);
    }
    if(isPointInGridBounds(bendBlock.x+1,bendBlock.y-1)){
        updateNetSurroundingNumber(2,netNumber,chipGridMatrix[bendBlock.x+1][bendBlock.y-1]);
    }
    if(isPointInGridBounds(bendBlock.x+1,bendBlock.y+1)){
        updateNetSurroundingNumber(2,netNumber,chipGridMatrix[bendBlock.x+1][bendBlock.y+1]);
    }
    if(isPointInGridBounds(bendBlock.x-1,bendBlock.y+1)){
        updateNetSurroundingNumber(2,netNumber,chipGridMatrix[bendBlock.x-1][bendBlock.y+1]);
    }
    if(isPointInGridBounds(bendBlock.x-1,bendBlock.y-1)){
        updateNetSurroundingNumber(2,netNumber,chipGridMatrix[bendBlock.x-1][bendBlock.y-1]);
    }
}

void unBlockCellsAroundTheBendForMetal2(int netNumber,chipUnitBlock bendBlock,chipUnitBlock source,chipUnitBlock target){
    clearNetSurroundingNumber(2,netNumber,chipGridMatrix[bendBlock.x][bendBlock.y]);
    if(isPointInGridBounds(bendBlock.x+1,bendBlock.y) && !isBlockInSurrounding(chipGridMatrix[bendBlock.x+1][bendBlock.y],source,target)){
        clearNetSurroundingNumber(2,netNumber,chipGridMatrix[bendBlock.x+1][bendBlock.y]);
    }
    if(isPointInGridBounds(bendBlock.x-1,bendBlock.y) && !isBlockInSurrounding(chipGridMatrix[bendBlock.x-1][bendBlock.y],source,target)){
        clearNetSurroundingNumber(2,netNumber,chipGridMatrix[bendBlock.x-1][bendBlock.y]);
    }
    if(isPointInGridBounds(bendBlock.x,bendBlock.y-1) && !isBlockInSurrounding(chipGridMatrix[bendBlock.x][bendBlock.y-1],source,target)){
        clearNetSurroundingNumber(2,netNumber,chipGridMatrix[bendBlock.x][bendBlock.y-1]);
    }
    if(isPointInGridBounds(bendBlock.x,bendBlock.y+1) && !isBlockInSurrounding(chipGridMatrix[bendBlock.x][bendBlock.y+1],source,target)){
        clearNetSurroundingNumber(2,netNumber,chipGridMatrix[bendBlock.x][bendBlock.y+1]);
    }
    if(isPointInGridBounds(bendBlock.x+1,bendBlock.y-1) && !isBlockInSurrounding(chipGridMatrix[bendBlock.x+1][bendBlock.y-1],source,target)){
        clearNetSurroundingNumber(2,netNumber,chipGridMatrix[bendBlock.x+1][bendBlock.y-1]);
    }
    if(isPointInGridBounds(bendBlock.x+1,bendBlock.y+1) && !isBlockInSurrounding(chipGridMatrix[bendBlock.x+1][bendBlock.y+1],source,target)){
        clearNetSurroundingNumber(2,netNumber,chipGridMatrix[bendBlock.x+1][bendBlock.y+1]);
    }
    if(isPointInGridBounds(bendBlock.x-1,bendBlock.y+1) && !isBlockInSurrounding(chipGridMatrix[bendBlock.x-1][bendBlock.y+1],source,target)){
        clearNetSurroundingNumber(2,netNumber,chipGridMatrix[bendBlock.x-1][bendBlock.y+1]);
    }
    if(isPointInGridBounds(bendBlock.x-1,bendBlock.y-1) && !isBlockInSurrounding(chipGridMatrix[bendBlock.x-1][bendBlock.y-1],source,target)){
        clearNetSurroundingNumber(2,netNumber,chipGridMatrix[bendBlock.x-1][bendBlock.y-1]);
    }
}

void insertViaAtGivenLocation(int netNumber,chipUnitBlock viaBlock){
    chipGridMatrix[viaBlock.x][viaBlock.y].via = 1;
    chipGridMatrix[viaBlock.x][viaBlock.y].m1NetNumber = netNumber;
    chipGridMatrix[viaBlock.x][viaBlock.y].m2NetNumber = netNumber;
    chipGridMatrix[viaBlock.x][viaBlock.y].m1Blocked = 1;//via
    chipGridMatrix[viaBlock.x][viaBlock.y].m2Blocked = 1;//via
    if(isPointInGridBounds(viaBlock.x+1,viaBlock.y)){
        chipGridMatrix[viaBlock.x+1][viaBlock.y].via = 2;//via beside block
        chipGridMatrix[viaBlock.x+1][viaBlock.y].m1Blocked = 1;//via
        chipGridMatrix[viaBlock.x+1][viaBlock.y].m2Blocked = 1;//via
        chipGridMatrix[viaBlock.x+1][viaBlock.y].m1NetSurrounding = netNumber;
        chipGridMatrix[viaBlock.x+1][viaBlock.y].m2NetSurrounding = netNumber;
    }
    if(isPointInGridBounds(viaBlock.x-1,viaBlock.y)){
        chipGridMatrix[viaBlock.x-1][viaBlock.y].via = 2;
        chipGridMatrix[viaBlock.x-1][viaBlock.y].m1Blocked = 1;//via
        chipGridMatrix[viaBlock.x-1][viaBlock.y].m2Blocked = 1;//via
        chipGridMatrix[viaBlock.x-1][viaBlock.y].m1NetSurrounding = netNumber;
        chipGridMatrix[viaBlock.x-1][viaBlock.y].m2NetSurrounding = netNumber;

    }
    if(isPointInGridBounds(viaBlock.x,viaBlock.y-1)){
        chipGridMatrix[viaBlock.x][viaBlock.y-1].via = 2;
        chipGridMatrix[viaBlock.x][viaBlock.y-1].m1Blocked = 1;//via
        chipGridMatrix[viaBlock.x][viaBlock.y-1].m2Blocked = 1;//via
        chipGridMatrix[viaBlock.x][viaBlock.y-1].m1NetSurrounding = netNumber;
        chipGridMatrix[viaBlock.x][viaBlock.y-1].m2NetSurrounding = netNumber;
    }
    if(isPointInGridBounds(viaBlock.x,viaBlock.y+1)){
        chipGridMatrix[viaBlock.x][viaBlock.y+1].via = 2;
        chipGridMatrix[viaBlock.x][viaBlock.y+1].m1Blocked = 1;//via
        chipGridMatrix[viaBlock.x][viaBlock.y+1].m2Blocked = 1;//via
        chipGridMatrix[viaBlock.x][viaBlock.y+1].m1NetSurrounding = netNumber;
        chipGridMatrix[viaBlock.x][viaBlock.y+1].m2NetSurrounding = netNumber;
    }
    if(isPointInGridBounds(viaBlock.x+1,viaBlock.y-1)){
        chipGridMatrix[viaBlock.x+1][viaBlock.y-1].via = 2;
        chipGridMatrix[viaBlock.x+1][viaBlock.y-1].m1Blocked = 1;//via
        chipGridMatrix[viaBlock.x+1][viaBlock.y-1].m2Blocked = 1;//via
        chipGridMatrix[viaBlock.x+1][viaBlock.y-1].m1NetSurrounding = netNumber;
        chipGridMatrix[viaBlock.x+1][viaBlock.y-1].m2NetSurrounding = netNumber;
    }
    if(isPointInGridBounds(viaBlock.x+1,viaBlock.y+1)){
        chipGridMatrix[viaBlock.x+1][viaBlock.y+1].via = 2;
        chipGridMatrix[viaBlock.x+1][viaBlock.y+1].m1Blocked = 1;//via
        chipGridMatrix[viaBlock.x+1][viaBlock.y+1].m2Blocked = 1;//via
        chipGridMatrix[viaBlock.x+1][viaBlock.y+1].m1NetSurrounding = netNumber;
        chipGridMatrix[viaBlock.x+1][viaBlock.y+1].m2NetSurrounding = netNumber;
    }
    if(isPointInGridBounds(viaBlock.x-1,viaBlock.y+1)){
        chipGridMatrix[viaBlock.x-1][viaBlock.y+1].via = 2;
        chipGridMatrix[viaBlock.x-1][viaBlock.y+1].m1Blocked = 1;//via
        chipGridMatrix[viaBlock.x-1][viaBlock.y+1].m2Blocked = 1;//via
        chipGridMatrix[viaBlock.x-1][viaBlock.y+1].m1NetSurrounding = netNumber;
        chipGridMatrix[viaBlock.x-1][viaBlock.y+1].m2NetSurrounding = netNumber;

    }
    if(isPointInGridBounds(viaBlock.x-1,viaBlock.y-1)){
        chipGridMatrix[viaBlock.x-1][viaBlock.y-1].via = 2;
        chipGridMatrix[viaBlock.x-1][viaBlock.y-1].m1Blocked = 1;//via
        chipGridMatrix[viaBlock.x-1][viaBlock.y-1].m2Blocked = 1;//via
        chipGridMatrix[viaBlock.x-1][viaBlock.y-1].m1NetSurrounding = netNumber;
        chipGridMatrix[viaBlock.x-1][viaBlock.y-1].m2NetSurrounding = netNumber;
    }
}

void deleteViaAtGivenLocation(chipUnitBlock viaBlock){
    chipGridMatrix[viaBlock.x][viaBlock.y].m1Blocked = 0;//via
    chipGridMatrix[viaBlock.x][viaBlock.y].m2Blocked = 0;//via
    chipGridMatrix[viaBlock.x][viaBlock.y].via = 0;
    if(isPointInGridBounds(viaBlock.x+1,viaBlock.y)){
        chipGridMatrix[viaBlock.x+1][viaBlock.y].via = 0;//via beside block
    }
    if(isPointInGridBounds(viaBlock.x-1,viaBlock.y)){
        chipGridMatrix[viaBlock.x-1][viaBlock.y].via = 0;

    }
    if(isPointInGridBounds(viaBlock.x,viaBlock.y-1)){
        chipGridMatrix[viaBlock.x][viaBlock.y-1].via = 0;
    }
    if(isPointInGridBounds(viaBlock.x,viaBlock.y+1)){
        chipGridMatrix[viaBlock.x][viaBlock.y+1].via = 0;
    }
    if(isPointInGridBounds(viaBlock.x+1,viaBlock.y-1)){
        chipGridMatrix[viaBlock.x+1][viaBlock.y-1].via = 0;
    }
    if(isPointInGridBounds(viaBlock.x+1,viaBlock.y+1)){
        chipGridMatrix[viaBlock.x+1][viaBlock.y+1].via = 0;
    }
    if(isPointInGridBounds(viaBlock.x-1,viaBlock.y+1)){
        chipGridMatrix[viaBlock.x-1][viaBlock.y+1].via = 0;
    }
    if(isPointInGridBounds(viaBlock.x-1,viaBlock.y-1)){
        chipGridMatrix[viaBlock.x-1][viaBlock.y-1].via = 0;
    }
}

bool checkIfViaCanBeInserted(chipUnitBlock viaBlock){
    bool viaInsert = true;
    if(isPointInGridBounds(viaBlock.x+1,viaBlock.y)){
        viaInsert = viaInsert && chipGridMatrix[viaBlock.x+1][viaBlock.y].via == 0;//via beside block
    }
    if(isPointInGridBounds(viaBlock.x-1,viaBlock.y)){
        viaInsert = viaInsert && chipGridMatrix[viaBlock.x-1][viaBlock.y].via == 0;

    }
    if(isPointInGridBounds(viaBlock.x,viaBlock.y-1)){
        viaInsert = viaInsert && chipGridMatrix[viaBlock.x][viaBlock.y-1].via  == 0;
    }
    if(isPointInGridBounds(viaBlock.x,viaBlock.y+1)){
        viaInsert = viaInsert && chipGridMatrix[viaBlock.x][viaBlock.y+1].via == 0;
    }
    if(isPointInGridBounds(viaBlock.x+1,viaBlock.y-1)){
        viaInsert = viaInsert && chipGridMatrix[viaBlock.x+1][viaBlock.y-1].via == 0;
    }
    if(isPointInGridBounds(viaBlock.x+1,viaBlock.y+1)){
        viaInsert = viaInsert && chipGridMatrix[viaBlock.x+1][viaBlock.y+1].via == 0;
    }
    if(isPointInGridBounds(viaBlock.x-1,viaBlock.y+1)){
        viaInsert = viaInsert && chipGridMatrix[viaBlock.x-1][viaBlock.y+1].via == 0;
    }
    if(isPointInGridBounds(viaBlock.x-1,viaBlock.y-1)){
        viaInsert = viaInsert && chipGridMatrix[viaBlock.x-1][viaBlock.y-1].via == 0;
    }
    return viaInsert;
}

int backTraceNetPath(int netNumber,chipUnitBlock source,chipUnitBlock target){
    int pathFound = -2;
    int targetLeeNumber = 0;
    chipGridMatrix[target.x][target.y].m1NetNumber = netNumber;
    int chipGridMatrixSize = chipGridMatrix.size();
    chipUnitBlock newTarget;
    chipUnitBlock m1FailedTarget;
    chipUnitBlock m2FailedTarget;
    unVisitedQueue.push(target);
    queue<chipUnitBlock> unVisitedM2Queue;
    if(target.m1LeeNumber== -1){
        targetLeeNumber = target.m2LeeNumber;
    }else{
        targetLeeNumber = target.m1LeeNumber;
    }
    //targetLeeNumber = target.leeNumber;
    bool popDone = false;
    bool fromWhileLoop = false;
    bool fromWhileLoopM2 = false;
    string previous = "";
    while(!unVisitedQueue.empty()){
        fromWhileLoop = false;
        popDone = false;
        newTarget = unVisitedQueue.front();
        if(newTarget.x == source.x && newTarget.y == source.y){
            pathFound = 0;
            break;
        }
        while(newTarget.y+1<chipGridMatrixSize && ((chipGridMatrix[newTarget.x][newTarget.y+1].m1Blocked == 0) || (chipGridMatrix[newTarget.x][newTarget.y+1].via == 2) || isBlockInSurrounding(chipGridMatrix[newTarget.x][newTarget.y+1],source,target) || (newTarget.x == source.x && newTarget.y+1 == source.y))
               && chipGridMatrix[newTarget.x][newTarget.y+1].m1LeeNumber == targetLeeNumber-1){
                chipGridMatrix[newTarget.x][newTarget.y+1].m1NetNumber = netNumber;
                updateNetSurroundingNumber(1,netNumber,chipGridMatrix[newTarget.x][newTarget.y+1]);
                unVisitedQueue.push(chipGridMatrix[newTarget.x][newTarget.y+1]);
                targetLeeNumber = targetLeeNumber-1;
                if(newTarget.x+1<chipGridMatrixSize){
                    updateNetSurroundingNumber(1,netNumber,chipGridMatrix[newTarget.x+1][newTarget.y+1]);
                }
                if(newTarget.x-1>=0){
                    updateNetSurroundingNumber(1,netNumber,chipGridMatrix[newTarget.x-1][newTarget.y+1]);
                }
                if(previous == "top" || previous == "bottom"){
                    blockCellsAroundTheBendForMetal1(netNumber,newTarget);
                }
                unVisitedQueue.pop();
                popDone = true;
                fromWhileLoop = true;
                newTarget = unVisitedQueue.front();
                if(fromWhileLoopM2 && !(newTarget.x == source.x && newTarget.y == source.y)){
                    insertViaAtGivenLocation(netNumber,newTarget);
                    fromWhileLoopM2 = false;
                }
                previous = "right";
        }
        if(popDone && !(newTarget.x == source.x && newTarget.y == source.y)){
            if(newTarget.x+1<chipGridMatrixSize){
                clearNetSurroundingNumber(1,netNumber,chipGridMatrix[newTarget.x+1][newTarget.y]);
            }
            if(newTarget.x-1>=0){
                clearNetSurroundingNumber(1,netNumber,chipGridMatrix[newTarget.x-1][newTarget.y]);
            }
            popDone = false;
        }
        while(newTarget.x+1<chipGridMatrixSize && ((chipGridMatrix[newTarget.x+1][newTarget.y].m1Blocked ==0) || (chipGridMatrix[newTarget.x+1][newTarget.y].via == 2)|| isBlockInSurrounding(chipGridMatrix[newTarget.x+1][newTarget.y],source,target) || (newTarget.x+1 == source.x && newTarget.y == source.y))
               && chipGridMatrix[newTarget.x+1][newTarget.y].m1LeeNumber == targetLeeNumber-1){
                chipGridMatrix[newTarget.x+1][newTarget.y].m1NetNumber = netNumber;
                updateNetSurroundingNumber(1,netNumber,chipGridMatrix[newTarget.x+1][newTarget.y]);
                unVisitedQueue.push(chipGridMatrix[newTarget.x+1][newTarget.y]);
                targetLeeNumber = targetLeeNumber-1;
                if(newTarget.y+1<chipGridMatrixSize){
                    updateNetSurroundingNumber(1,netNumber,chipGridMatrix[newTarget.x+1][newTarget.y+1]);
                }
                if(newTarget.y-1>=0){
                    updateNetSurroundingNumber(1,netNumber,chipGridMatrix[newTarget.x+1][newTarget.y-1]);
                }
                if(previous == "left" || previous == "right"){
                    blockCellsAroundTheBendForMetal1(netNumber,newTarget);
                }
                unVisitedQueue.pop();
                popDone = true;
                fromWhileLoop = true;
                newTarget = unVisitedQueue.front();
                if(fromWhileLoopM2 && !(newTarget.x == source.x && newTarget.y == source.y)){
                    insertViaAtGivenLocation(netNumber,newTarget);
                    fromWhileLoopM2 = false;
                }
                previous = "bottom";
        }
        if(popDone && !(newTarget.x == source.x && newTarget.y == source.y)){
            if(newTarget.y+1<chipGridMatrixSize){
                clearNetSurroundingNumber(1,netNumber,chipGridMatrix[newTarget.x][newTarget.y+1]);
            }
            if(newTarget.y-1>=0){
                clearNetSurroundingNumber(1,netNumber,chipGridMatrix[newTarget.x][newTarget.y-1]);
            }
            popDone = false;
        }
        while(newTarget.x-1>=0 && ((chipGridMatrix[newTarget.x-1][newTarget.y].m1Blocked==0) || (chipGridMatrix[newTarget.x-1][newTarget.y].via==2)|| isBlockInSurrounding(chipGridMatrix[newTarget.x-1][newTarget.y],source,target) || (newTarget.x-1 == source.x && newTarget.y == source.y))
               && chipGridMatrix[newTarget.x-1][newTarget.y].m1LeeNumber == targetLeeNumber-1){
                chipGridMatrix[newTarget.x-1][newTarget.y].m1NetNumber = netNumber;
                updateNetSurroundingNumber(1,netNumber,chipGridMatrix[newTarget.x-1][newTarget.y]);
                unVisitedQueue.push(chipGridMatrix[newTarget.x-1][newTarget.y]);
                targetLeeNumber = targetLeeNumber-1;
                if(newTarget.y+1<chipGridMatrixSize){
                    updateNetSurroundingNumber(1,netNumber,chipGridMatrix[newTarget.x-1][newTarget.y+1]);
                }
                if(newTarget.y-1>=0){
                    updateNetSurroundingNumber(1,netNumber,chipGridMatrix[newTarget.x-1][newTarget.y-1]);
                }
                if(previous == "left" || previous == "right"){
                    blockCellsAroundTheBendForMetal1(netNumber,newTarget);
                }
                /*if(newTarget.x+1<chipGridMatrixSize){
                    chipGridMatrix[newTarget.x+1][newTarget.y].m1Blocked = 1;
                }*/
                unVisitedQueue.pop();
                popDone = true;
                fromWhileLoop = true;
                newTarget = unVisitedQueue.front();
                if(fromWhileLoopM2 && !(newTarget.x == source.x && newTarget.y == source.y)){
                    insertViaAtGivenLocation(netNumber,newTarget);
                    fromWhileLoopM2 = false;
                }
                previous = "top";
        }
        if(popDone && !(newTarget.x == source.x && newTarget.y == source.y)){
            if(newTarget.y+1<chipGridMatrixSize){
                clearNetSurroundingNumber(1,netNumber,chipGridMatrix[newTarget.x][newTarget.y+1]);
            }
            if(newTarget.y-1>=0){
                clearNetSurroundingNumber(1,netNumber,chipGridMatrix[newTarget.x][newTarget.y-1]);
            }
            popDone = false;
        }
        while(newTarget.y-1>=0 && ((chipGridMatrix[newTarget.x][newTarget.y-1].m1Blocked == 0) || (chipGridMatrix[newTarget.x][newTarget.y-1].via == 2) || isBlockInSurrounding(chipGridMatrix[newTarget.x][newTarget.y-1],source,target) || (newTarget.x == source.x && newTarget.y-1 == source.y))
               && chipGridMatrix[newTarget.x][newTarget.y-1].m1LeeNumber == targetLeeNumber-1){
                chipGridMatrix[newTarget.x][newTarget.y-1].m1NetNumber = netNumber;
                updateNetSurroundingNumber(1,netNumber,chipGridMatrix[newTarget.x][newTarget.y-1]);
                unVisitedQueue.push(chipGridMatrix[newTarget.x][newTarget.y-1]);
                targetLeeNumber = targetLeeNumber-1;
                if(newTarget.x+1<chipGridMatrixSize){
                    updateNetSurroundingNumber(1,netNumber,chipGridMatrix[newTarget.x+1][newTarget.y-1]);
                }
                if(newTarget.x-1>=0){
                    updateNetSurroundingNumber(1,netNumber,chipGridMatrix[newTarget.x-1][newTarget.y-1]);
                }
                if(previous == "top" || previous == "bottom"){
                    blockCellsAroundTheBendForMetal1(netNumber,newTarget);
                }
                /*if(newTarget.y+1<chipGridMatrixSize){
                    chipGridMatrix[newTarget.x][newTarget.y+1].m1Blocked = 1;
                }*/
                unVisitedQueue.pop();
                popDone = true;
                fromWhileLoop = true;
                newTarget = unVisitedQueue.front();
                if(fromWhileLoopM2 && !(newTarget.x == source.x && newTarget.y == source.y)){
                    insertViaAtGivenLocation(netNumber,newTarget);
                    fromWhileLoopM2 = false;
                }
                previous = "left";
        }
        if(popDone && !(newTarget.x == source.x && newTarget.y == source.y)){
            if(newTarget.x+1<chipGridMatrixSize){
                clearNetSurroundingNumber(1,netNumber,chipGridMatrix[newTarget.x+1][newTarget.y]);
            }
            if(newTarget.x-1>=0){
                clearNetSurroundingNumber(1,netNumber,chipGridMatrix[newTarget.x-1][newTarget.y]);
            }
            popDone = false;
        }
        if(!fromWhileLoop){
            m1FailedTarget = newTarget;

            /*if(previous == "right"){
                insertViaAtGivenLocation(netNumber,chipGridMatrix[newTarget.x][newTarget.y-1]);
            }else if(previous == "left"){
                insertViaAtGivenLocation(netNumber,chipGridMatrix[newTarget.x][newTarget.y+1]);
            }else if(previous == "top"){
                insertViaAtGivenLocation(netNumber,chipGridMatrix[newTarget.x+1][newTarget.y]);
            }else if(previous == "bottom"){
                insertViaAtGivenLocation(netNumber,chipGridMatrix[newTarget.x-1][newTarget.y]);
            }*/
            if(!(newTarget.x == source.x && newTarget.y == source.y)){
                insertViaAtGivenLocation(netNumber,newTarget);
            }
            popDone = false;
            fromWhileLoopM2 = false;
            unVisitedM2Queue.push(newTarget);
            while(!unVisitedM2Queue.empty()){
                fromWhileLoopM2 = false;
                popDone = false;
            //if(previous == "right"){
                //printChipGridMatrixElement(newTarget.x,newTarget.y+1);
                while(newTarget.y+1<chipGridMatrixSize && ((chipGridMatrix[newTarget.x][newTarget.y+1].m2Blocked == 0) || (chipGridMatrix[newTarget.x][newTarget.y+1].via == 2) || isBlockInSurrounding(chipGridMatrix[newTarget.x][newTarget.y+1],source,target) || (newTarget.x == source.x && newTarget.y+1 == source.y))
                   && chipGridMatrix[newTarget.x][newTarget.y+1].m1LeeNumber == -1 && chipGridMatrix[newTarget.x][newTarget.y+1].m2LeeNumber == targetLeeNumber-1){
                    chipGridMatrix[newTarget.x][newTarget.y+1].m2NetNumber = netNumber;
                    updateNetSurroundingNumber(2,netNumber,chipGridMatrix[newTarget.x][newTarget.y+1]);
                    unVisitedQueue.push(chipGridMatrix[newTarget.x][newTarget.y+1]);
                    unVisitedM2Queue.push(chipGridMatrix[newTarget.x][newTarget.y+1]);
                    targetLeeNumber = targetLeeNumber-1;
                    if(newTarget.x+1<chipGridMatrixSize){
                        updateNetSurroundingNumber(2,netNumber,chipGridMatrix[newTarget.x+1][newTarget.y+1]);
                    }
                    if(newTarget.x-1>=0){
                        updateNetSurroundingNumber(2,netNumber,chipGridMatrix[newTarget.x-1][newTarget.y+1]);
                    }
                    if(previous == "top" || previous == "bottom"){
                        blockCellsAroundTheBendForMetal2(netNumber,newTarget);
                    }
                    unVisitedQueue.pop();
                    unVisitedM2Queue.pop();
                    popDone = true;
                    //fromWhileLoop = true;
                    fromWhileLoopM2 = true;
                    newTarget = unVisitedQueue.front();
                    previous = "right";
                }
                if(popDone && !(newTarget.x == source.x && newTarget.y == source.y)){
                    if(newTarget.x+1<chipGridMatrixSize){
                        clearNetSurroundingNumber(2,netNumber,chipGridMatrix[newTarget.x+1][newTarget.y]);
                    }
                    if(newTarget.x-1>=0){
                        clearNetSurroundingNumber(2,netNumber,chipGridMatrix[newTarget.x-1][newTarget.y]);
                    }
                    popDone = false;
                }
            //}else if(previous == "left"){
                //printChipGridMatrixElement(newTarget.x,newTarget.y-1);
                while(newTarget.y-1>=0 && ((chipGridMatrix[newTarget.x][newTarget.y-1].m2Blocked == 0) || (chipGridMatrix[newTarget.x][newTarget.y-1].via == 2) || isBlockInSurrounding(chipGridMatrix[newTarget.x][newTarget.y-1],source,target) || (newTarget.x == source.x && newTarget.y-1 == source.y))
                   && chipGridMatrix[newTarget.x][newTarget.y-1].m1LeeNumber == -1 && chipGridMatrix[newTarget.x][newTarget.y-1].m2LeeNumber == targetLeeNumber-1){
                    chipGridMatrix[newTarget.x][newTarget.y-1].m2NetNumber = netNumber;
                    updateNetSurroundingNumber(2,netNumber,chipGridMatrix[newTarget.x][newTarget.y-1]);
                    unVisitedQueue.push(chipGridMatrix[newTarget.x][newTarget.y-1]);
                    unVisitedM2Queue.push(chipGridMatrix[newTarget.x][newTarget.y-1]);
                    targetLeeNumber = targetLeeNumber-1;
                    if(newTarget.x+1<chipGridMatrixSize){
                        updateNetSurroundingNumber(2,netNumber,chipGridMatrix[newTarget.x+1][newTarget.y-1]);
                    }
                    if(newTarget.x-1>=0){
                        updateNetSurroundingNumber(2,netNumber,chipGridMatrix[newTarget.x-1][newTarget.y-1]);
                    }
                    if(previous == "top" || previous == "bottom"){
                        blockCellsAroundTheBendForMetal2(netNumber,newTarget);
                    }
                    /*if(newTarget.y+1<chipGridMatrixSize){
                        chipGridMatrix[newTarget.x][newTarget.y+1].m1Blocked = 1;
                    }*/
                    unVisitedQueue.pop();
                    unVisitedM2Queue.pop();
                    popDone = true;
                    //fromWhileLoop = true;
                    fromWhileLoopM2 = true;
                    newTarget = unVisitedQueue.front();
                    previous = "left";
                }
                if(popDone && !(newTarget.x == source.x && newTarget.y == source.y)){
                    if(newTarget.x+1<chipGridMatrixSize){
                        clearNetSurroundingNumber(2,netNumber,chipGridMatrix[newTarget.x+1][newTarget.y]);
                    }
                    if(newTarget.x-1>=0){
                        clearNetSurroundingNumber(2,netNumber,chipGridMatrix[newTarget.x-1][newTarget.y]);
                    }
                    popDone = false;
                }
            //}else if(previous == "top"){
                //printChipGridMatrixElement(newTarget.x-1,newTarget.y-1);
                while(newTarget.x-1>=0 && ((chipGridMatrix[newTarget.x-1][newTarget.y].m2Blocked == 0) || (chipGridMatrix[newTarget.x-1][newTarget.y].via == 2) || isBlockInSurrounding(chipGridMatrix[newTarget.x-1][newTarget.y],source,target) || (newTarget.x-1 == source.x && newTarget.y == source.y))
                   && chipGridMatrix[newTarget.x-1][newTarget.y].m1LeeNumber == -1 && chipGridMatrix[newTarget.x-1][newTarget.y].m2LeeNumber == targetLeeNumber-1){
                    chipGridMatrix[newTarget.x-1][newTarget.y].m2NetNumber = netNumber;
                    updateNetSurroundingNumber(2,netNumber,chipGridMatrix[newTarget.x-1][newTarget.y]);
                    unVisitedQueue.push(chipGridMatrix[newTarget.x-1][newTarget.y]);
                    unVisitedM2Queue.push(chipGridMatrix[newTarget.x-1][newTarget.y]);
                    targetLeeNumber = targetLeeNumber-1;
                    if(newTarget.y+1<chipGridMatrixSize){
                        updateNetSurroundingNumber(2,netNumber,chipGridMatrix[newTarget.x-1][newTarget.y+1]);
                    }
                    if(newTarget.y-1>=0){
                        updateNetSurroundingNumber(2,netNumber,chipGridMatrix[newTarget.x-1][newTarget.y-1]);
                    }
                    if(previous == "left" || previous == "right"){
                        blockCellsAroundTheBendForMetal2(netNumber,newTarget);
                    }
                    /*if(newTarget.x+1<chipGridMatrixSize){
                        chipGridMatrix[newTarget.x+1][newTarget.y].m1Blocked = 1;
                    }*/
                    unVisitedQueue.pop();
                    unVisitedM2Queue.pop();
                    popDone = true;
                    //fromWhileLoop = true;
                    fromWhileLoopM2 = true;
                    newTarget = unVisitedQueue.front();
                    previous = "top";
                }
                if(popDone && !(newTarget.x == source.x && newTarget.y == source.y)){
                    if(newTarget.y+1<chipGridMatrixSize){
                        clearNetSurroundingNumber(2,netNumber,chipGridMatrix[newTarget.x][newTarget.y+1]);
                    }
                    if(newTarget.y-1>=0){
                        clearNetSurroundingNumber(2,netNumber,chipGridMatrix[newTarget.x][newTarget.y-1]);
                    }
                    popDone = false;
                }
            //}else if(previous == "bottom"){
                //printChipGridMatrixElement(newTarget.x+1,newTarget.y);
                while(newTarget.x+1<chipGridMatrixSize && ((chipGridMatrix[newTarget.x+1][newTarget.y].m2Blocked ==0) || (chipGridMatrix[newTarget.x+1][newTarget.y].via == 2) || isBlockInSurrounding(chipGridMatrix[newTarget.x+1][newTarget.y],source,target) || (newTarget.x+1 == source.x && newTarget.y == source.y))
                   && chipGridMatrix[newTarget.x+1][newTarget.y].m1LeeNumber == -1 && chipGridMatrix[newTarget.x+1][newTarget.y].m2LeeNumber == targetLeeNumber-1){
                    chipGridMatrix[newTarget.x+1][newTarget.y].m2NetNumber = netNumber;
                    updateNetSurroundingNumber(2,netNumber,chipGridMatrix[newTarget.x+1][newTarget.y]);
                    unVisitedQueue.push(chipGridMatrix[newTarget.x+1][newTarget.y]);
                    unVisitedM2Queue.push(chipGridMatrix[newTarget.x+1][newTarget.y]);
                    targetLeeNumber = targetLeeNumber-1;
                    if(newTarget.y+1<chipGridMatrixSize){
                        updateNetSurroundingNumber(2,netNumber,chipGridMatrix[newTarget.x+1][newTarget.y+1]);
                    }
                    if(newTarget.y-1>=0){
                        updateNetSurroundingNumber(2,netNumber,chipGridMatrix[newTarget.x+1][newTarget.y-1]);
                    }
                    if(previous == "left" || previous == "right"){
                        blockCellsAroundTheBendForMetal2(netNumber,newTarget);
                    }
                    unVisitedQueue.pop();
                    unVisitedM2Queue.pop();
                    popDone = true;
                    //fromWhileLoop = true;
                    fromWhileLoopM2 = true;
                    newTarget = unVisitedQueue.front();
                    previous = "bottom";
                }
                if(popDone && !(newTarget.x == source.x && newTarget.y == source.y)){
                    if(newTarget.y+1<chipGridMatrixSize){
                        clearNetSurroundingNumber(2,netNumber,chipGridMatrix[newTarget.x][newTarget.y+1]);
                    }
                    if(newTarget.y-1>=0){
                        clearNetSurroundingNumber(2,netNumber,chipGridMatrix[newTarget.x][newTarget.y-1]);
                    }
                    popDone = false;
                }
            //}
            //unVisitedQueue.pop();
                if(!fromWhileLoopM2){
                    m2FailedTarget = newTarget;
                    fromWhileLoopM2 = true;
                    unVisitedM2Queue.pop();
                }
                //fromWhileLoop = true;
            }

        }
        if(m1FailedTarget.x == m2FailedTarget.x && m1FailedTarget.y == m2FailedTarget.y){
            unVisitedQueue.pop();
        }
    }
    chipGridMatrix[source.x][source.y].m1NetNumber = netNumber;
    while(!unVisitedQueue.empty()){
        unVisitedQueue.pop();
    }
    return pathFound;
}

void clearLeeNumbersFromSourceToTerminal(int netNumber,chipUnitBlock source,chipUnitBlock target){
    int sourceX = source.x;
    int sourceY = source.y;
    int chipGridMatrixSize = chipGridMatrix.size();
    int currentX = -1;
    int currentY = -1;
    chipUnitBlock newSource;
    unVisitedQueue.push(source);
    while(!unVisitedQueue.empty()){
        newSource = unVisitedQueue.front();
        sourceX = newSource.x;
        sourceY = newSource.y;
        //chipGridMatrix[sourceX][sourceY].leeNumber = -1;
        chipGridMatrix[sourceX][sourceY].m1LeeNumber = -1;
        chipGridMatrix[sourceX][sourceY].m2LeeNumber = -1;
        if(newSource.x == target.x && newSource.y == target.y){
            break;
        }
        currentX = sourceX+1;
        currentY = sourceY;
        if(currentX<chipGridMatrixSize && (chipGridMatrix[currentX][currentY].m1LeeNumber != -1 || chipGridMatrix[currentX][currentY].m2LeeNumber != -1)){
            //chipGridMatrix[currentX][currentY].leeNumber = -1;
            chipGridMatrix[currentX][currentY].m1LeeNumber = -1;
            chipGridMatrix[currentX][currentY].m2LeeNumber = -1;
            unVisitedQueue.push(chipGridMatrix[currentX][currentY]);
        }
        currentX = sourceX-1;
        currentY = sourceY;
        if(currentX>=0 && (chipGridMatrix[currentX][currentY].m1LeeNumber != -1 || chipGridMatrix[currentX][currentY].m2LeeNumber != -1)){
            //chipGridMatrix[currentX][currentY].leeNumber = -1;
            chipGridMatrix[currentX][currentY].m1LeeNumber = -1;
            chipGridMatrix[currentX][currentY].m2LeeNumber = -1;
            unVisitedQueue.push(chipGridMatrix[currentX][currentY]);
        }
        currentX = sourceX;
        currentY = sourceY+1;
        if(currentY<chipGridMatrixSize && (chipGridMatrix[currentX][currentY].m1LeeNumber != -1 || chipGridMatrix[currentX][currentY].m2LeeNumber != -1)){
            //chipGridMatrix[currentX][currentY].leeNumber = -1;
            chipGridMatrix[currentX][currentY].m1LeeNumber = -1;
            chipGridMatrix[currentX][currentY].m2LeeNumber = -1;
            unVisitedQueue.push(chipGridMatrix[currentX][currentY]);
        }
        currentX = sourceX;
        currentY = sourceY-1;
        if(currentY>=0 && (chipGridMatrix[currentX][currentY].m1LeeNumber != -1 || chipGridMatrix[currentX][currentY].m2LeeNumber != -1)){
            //chipGridMatrix[currentX][currentY].leeNumber = -1;
            chipGridMatrix[currentX][currentY].m1LeeNumber = -1;
            chipGridMatrix[currentX][currentY].m2LeeNumber = -1;
            unVisitedQueue.push(chipGridMatrix[currentX][currentY]);
        }
        unVisitedQueue.pop();
    }
    while(!unVisitedQueue.empty()){
        newSource = unVisitedQueue.front();
        sourceX = newSource.x;
        sourceY = newSource.y;
        //chipGridMatrix[sourceX][sourceY].leeNumber = -1;
        chipGridMatrix[sourceX][sourceY].m1LeeNumber = -1;
        chipGridMatrix[sourceX][sourceY].m2LeeNumber = -1;
        unVisitedQueue.pop();
    }
}
void convertChipGridMatrixToMagicFile(ofstream &outputMagicFile){
    int chipGridMatrixRows = chipGridMatrix.size();
    int chipGridMatrixColumns = -1;
    outputMagicFile<<"magic"<< endl;
    outputMagicFile<<"tech scmos"<< endl;
    outputMagicFile<<"timestamp 11111111111"<< endl;
    cellcoOrdinates cellPoint;
    int noOfCells = 0;
    int wireLength = 0;
    int chipArea = 0;
    int viaCount = 0;
    sort(netTerminalVector.begin(),netTerminalVector.end(),sortByNetNumber);
    for(int i=0;i<chipGridMatrixRows;i++){
        chipGridMatrixColumns = chipGridMatrix[i].size();
        for(int j=0;j<chipGridMatrixColumns;j++){
            //chipArea++;
            if(chipGridMatrix[i][j].m1NetNumber != 0){
                outputMagicFile<<"<< m1 >>"<<endl;
                outputMagicFile<<"rect "<<chipGridMatrix[i][j].y <<" "<<chipGridMatrix[i][j].x<<" "<<chipGridMatrix[i][j].y+1<<" "<<chipGridMatrix[i][j].x+1<<" "<<endl;
                wireLength++;
            }
            if(chipGridMatrix[i][j].m2NetNumber != 0){
                outputMagicFile<<"<< m2 >>"<<endl;
                outputMagicFile<<"rect "<<chipGridMatrix[i][j].y <<" "<<chipGridMatrix[i][j].x<<" "<<chipGridMatrix[i][j].y+1<<" "<<chipGridMatrix[i][j].x+1<<" "<<endl;
                wireLength++;
            }
            if(chipGridMatrix[i][j].via == 1 && (chipGridMatrix[i][j].m1NetNumber != 0  && chipGridMatrix[i][j].m2NetNumber != 0)){
                outputMagicFile<<"<< m2c >>"<<endl;
                outputMagicFile<<"rect "<<chipGridMatrix[i][j].y <<" "<<chipGridMatrix[i][j].x<<" "<<chipGridMatrix[i][j].y+1<<" "<<chipGridMatrix[i][j].x+1<<" "<<endl;
                outputMagicFile<<"<< m1 >>"<<endl;
                outputMagicFile<<"rect "<<chipGridMatrix[i][j].y <<" "<<chipGridMatrix[i][j].x<<" "<<chipGridMatrix[i][j].y+1<<" "<<chipGridMatrix[i][j].x+1<<" "<<endl;
                outputMagicFile<<"<< m2 >>"<<endl;
                outputMagicFile<<"rect "<<chipGridMatrix[i][j].y <<" "<<chipGridMatrix[i][j].x<<" "<<chipGridMatrix[i][j].y+1<<" "<<chipGridMatrix[i][j].x+1<<" "<<endl;
                wireLength++;
                viaCount++;
            }
            if(chipGridMatrix[i][j].cellNumber != -1){
                outputMagicFile<<"<< pdiffusion >>"<<endl;
                outputMagicFile<<"rect "<<chipGridMatrix[i][j].y<<" "<<chipGridMatrix[i][j].x<<" "<<chipGridMatrix[i][j].y+1<<" "<<chipGridMatrix[i][j].x+1<<" "<<endl;
            }
            /*if(chipGridMatrix[i][j].cellNumber != -1){
                outputMagicFile<<"<< labels >>"<<endl;
                outputMagicFile<<"rlabel "<<"pdiffusion "<<chipGridMatrix[i][j].x<<" "<<chipGridMatrix[i][j].y<<" "<<chipGridMatrix[i][j].x+1<<" "<<chipGridMatrix[i][j].y+1<<" "<<"0 "<<"cell no = "<<chipGridMatrix[i][j].cellNumber+1<<endl;
            }*/

        }
    }
    noOfCells = cellCoOrdinateVector.size();
    for(int i=0;i<noOfCells;i++){
        cellPoint = cellCoOrdinateVector.at(i);
        //if(cellPlacementMatrix[i][j].cellNumber != -1){
        outputMagicFile<<"<< labels >>"<<endl;
        outputMagicFile<<"rlabel "<<"pdiffusion "<<cellPoint.t1.y<<" "<<cellPoint.t1.x<<" "<<cellPoint.t1.y+1<<" "<<cellPoint.t1.x+1<<" "<<" 0 "<<"t = 1"<<endl;
        outputMagicFile<<"rlabel "<<"pdiffusion "<<cellPoint.t2.y<<" "<<cellPoint.t2.x<<" "<<cellPoint.t2.y+1<<" "<<cellPoint.t2.x+1<<" "<<" 0 "<<"t = 2"<<endl;
        outputMagicFile<<"rlabel "<<"pdiffusion "<<cellPoint.t3.y<< " " <<cellPoint.t3.x<< " " <<cellPoint.t3.y+1<<" "<<cellPoint.t3.x+1<<" "<<" 0 "<<"t = 3"<<endl;
        outputMagicFile<<"rlabel "<<"pdiffusion "<<cellPoint.t4.y<< " " <<cellPoint.t4.x<< " " <<cellPoint.t4.y+1<< " " <<cellPoint.t4.x+1<<" "<<" 0 "<<"t = 4"<<endl;
        outputMagicFile<<"rlabel "<<"pdiffusion "<<cellPoint.referencePoint.y<<" "<<cellPoint.referencePoint.x<<" "<<cellPoint.referencePoint.y+6<<" "<<cellPoint.referencePoint.x+6<<" "<<"0 "<<"cell no = "<<i+1<<endl;
        //}
        outputMagicFile<<"<< m1 >>"<<endl;
        outputMagicFile<<"rect "<<cellPoint.t1.y <<" "<<cellPoint.t1.x<<" "<<cellPoint.t1.y+1<<" "<<cellPoint.t1.x+1<<" "<<endl;
        outputMagicFile<<"rect "<<cellPoint.t2.y <<" "<<cellPoint.t2.x<<" "<<cellPoint.t2.y+1<<" "<<cellPoint.t2.x+1<<" "<<endl;
        outputMagicFile<<"rect "<<cellPoint.t3.y <<" "<<cellPoint.t3.x<<" "<<cellPoint.t3.y+1<<" "<<cellPoint.t3.x+1<<" "<<endl;
        outputMagicFile<<"rect "<<cellPoint.t4.y <<" "<<cellPoint.t4.x<<" "<<cellPoint.t4.y+1<<" "<<cellPoint.t4.x+1<<" "<<endl;

        outputMagicFile<<"<< m2 >>"<<endl;
        outputMagicFile<<"rect "<<cellPoint.t1.y <<" "<<cellPoint.t1.x<<" "<<cellPoint.t1.y+1<<" "<<cellPoint.t1.x+1<<" "<<endl;
        outputMagicFile<<"rect "<<cellPoint.t2.y <<" "<<cellPoint.t2.x<<" "<<cellPoint.t2.y+1<<" "<<cellPoint.t2.x+1<<" "<<endl;
        outputMagicFile<<"rect "<<cellPoint.t3.y <<" "<<cellPoint.t3.x<<" "<<cellPoint.t3.y+1<<" "<<cellPoint.t3.x+1<<" "<<endl;
        outputMagicFile<<"rect "<<cellPoint.t4.y <<" "<<cellPoint.t4.x<<" "<<cellPoint.t4.y+1<<" "<<cellPoint.t4.x+1<<" "<<endl;

        outputMagicFile<<"<< m2c >>"<<endl;
        outputMagicFile<<"rect "<<cellPoint.t1.y <<" "<<cellPoint.t1.x<<" "<<cellPoint.t1.y+1<<" "<<cellPoint.t1.x+1<<" "<<endl;
        outputMagicFile<<"rect "<<cellPoint.t2.y <<" "<<cellPoint.t2.x<<" "<<cellPoint.t2.y+1<<" "<<cellPoint.t2.x+1<<" "<<endl;
        outputMagicFile<<"rect "<<cellPoint.t3.y <<" "<<cellPoint.t3.x<<" "<<cellPoint.t3.y+1<<" "<<cellPoint.t3.x+1<<" "<<endl;
        outputMagicFile<<"rect "<<cellPoint.t4.y <<" "<<cellPoint.t4.x<<" "<<cellPoint.t4.y+1<<" "<<cellPoint.t4.x+1<<" "<<endl;
        viaCount++;
    }
    outputMagicFile<<"<< end >> "<<endl;
    outputMagicFile.close();
    cout<<"Final Wire length for routed nets: "<<wireLength<<endl;
    cout<<"No.of Vias: "<<viaCount<<endl;
    //cout<<"Chip Area : "<<chipArea<<endl;
}


