import sys
from matplotlib import pyplot as plt

"""
Point: defines a 2d point
    Parameter(s)
        x:  an int or a float represents x-coordinate value
        y:  an int or a float represents y-coordinate value
        flag:   indicate type of point it is
            P:  left endpoint
            Q:  right endpoint
"""
class Point:
    def __init__(self,x,y,flag):
        self.x=x
        self.y=y
        self.flag=flag

    def __eq__(self,p2):
        return self.x==p2.x and self.y==p2.y

    def __lt__(self,p2):
        if self.x<p2.x:
            return True
        elif self.x==p2.x:
            return self.y<p2.y
        return False
    
    def __gt__(self,p2):
        if self.x>p2.x:
            return True
        elif self.x==p2.x:
            return self.y>p2.y
        return False

    def __hash__(self):
        return hash(self.x+self.y)

    def __str__(self):
        return "(" + str(self.x) + ", " + str(self.y) + ")"

"""
LineSegment: Defines a line segment
    Parameter(s)
        start:  left endpoint of this line segment
        end:    right endpoint of this line segment
        m:  slope of this line segment equation
        b:  b of this line segment equation
"""
class LineSegment:
    def __init__(self,leftEndPoint,rightEndPoint):
        self.start=leftEndPoint
        self.end=rightEndPoint
        self.m=(leftEndPoint.y-rightEndPoint.y)/(leftEndPoint.x-rightEndPoint.x)
        self.b=leftEndPoint.y-self.m*leftEndPoint.x
    
    def __eq__(self,ls2):
        return self.start==ls2.start and self.end==ls2.end

    def get_y(self, x):
        return self.m * x + self.b

    def __str__(self):
        return self.start.__str__() + "-" + self.end.__str__()

"""
Wall: Defines a vertical line represents the wall created by endpoints
    Parameter(s)
        x:  x-coordinate value of this vertical line
        ys: y-coordinate value indicates the position of the point
        yi: the minimum y value this wall reaches
        ya: the maximum y value this wall reaches
"""
class Wall:
    def __init__(self,x,yStart,yMin,yMax):
        self.x=x
        self.ys=yStart
        self.yi=yMin
        self.ya=yMax
    
    def __eq__(self,w2):
        return self.x==w2.x
    
    def __lt__(self,w2):
        return self.x<w2.x
    
    def __gt__(self,w2):
        return self.x>w2.x

"""
Cell: Defines the trapezoid
    Parameter(s)
        left:   left boundary, can be either a Wall or a number(means the bounding box edge)
        right:  right boundary, can be either a Wall or a number(means the bounding box edge)
        above:  above boundary, can be either a line segment or a number(means the bounding box edge)
        below:  below boundary, can be either a line segment or a number(means the bounding box edge)
"""
class Cell:
    def __init__(self,leftBoundary,rightBoundary,aboveBoundary,belowBoundary):
        self.left=leftBoundary
        self.right=rightBoundary
        self.above=aboveBoundary
        self.below=belowBoundary
    
    def __eq__(self,c2):
        return self.left==c2.left and self.right==c2.right and self.above==c2.above and self.below==c2.below

"""
DAG: Defines a directed acyclic graph
    Parameter(s)
        leftChild:  left child of current node (or above child for LineSegment; None for Cell)
        itself: the actual object in this node, can be Point, LineSegment, or Cell
        rightChild: right child of current node (or below child for LineSegment; None for Cell)
"""
class DAG:
    def __init__(self,leftChild,current,rightChild):
        self.leftChild=leftChild
        self.itself=current
        self.rightChild=rightChild

    def printRecursive(self, indent):
        res = ""
        for i in range(indent):
            res += "\t"
        if isinstance(self.itself, Point):
            res += "Point " + self.itself.__str__()
        elif isinstance(self.itself, Cell):
            leftx = self.itself.left
            if isinstance(leftx, Wall):
                leftx = self.itself.left.x
            rightx = self.itself.right
            if isinstance(rightx, Wall):
                rightx = self.itself.right.x
            res += "Cell [x: " + str(leftx) + ", " + str(rightx) + "] [y: " + \
                self.itself.above.__str__() + ", " + self.itself.below.__str__() + "]"
        elif isinstance(self.itself, LineSegment):
            res += "Segment " + self.itself.__str__()
        print(res)
        if self.leftChild:
            self.leftChild.printRecursive(indent + 1)
        if self.rightChild:
            self.rightChild.printRecursive(indent + 1)

    def plot(self, filename):
        self.plotRecursive()
        plt.savefig(filename)

    def plotRecursive(self):
        if isinstance(self.itself, Point):
            plt.scatter([self.itself.x], [self.itself.y], c=[(0, 0, 0)])
        elif isinstance(self.itself, Cell):
            if(isinstance(self.itself.left, Wall)):
                plt.plot([self.itself.left.x, self.itself.left.x], [self.itself.left.yi, self.itself.left.ya], c=(0.5, 0.5, 1))
            if (isinstance(self.itself.right, Wall)):
                plt.plot([self.itself.right.x, self.itself.right.x], [self.itself.right.yi, self.itself.right.ya], c=(0.5, 0.5, 1))

        elif isinstance(self.itself, LineSegment):
            plt.plot([self.itself.start.x, self.itself.end.x], [self.itself.start.y, self.itself.end.y], c=(0, 0, 0))
        if self.leftChild:
            self.leftChild.plotRecursive()
        if self.rightChild:
            self.rightChild.plotRecursive()
"""
pointLocate: finds the trapezoid a point is located
    Parameter(s):
        p:  point
        node: current DAG node
"""
def pointLocate(p,node):
    #if current node is a Cell, return
    if isinstance(node.itself,Cell):
        return node
    #else if current node is a Point, search through its left or right child
    elif isinstance(node.itself,Point):
        if p<node.itself:
            return pointLocate(p,node.leftChild)
        else:
            return pointLocate(p,node.rightChild)
    #else if current node is a LineSegment, search through its above or below child
    elif isinstance(node.itself,LineSegment):
        y=node.itself.m*p.x+node.itself.b
        if p.y<y:
            return pointLocate(p,node.rightChild)
        else:
            return pointLocate(p,node.leftChild)

"""
wallCreate: creates walls generated by the endpoints of a line segment
"""
def wallCreate(pCell,qCell,p,q):
    #check left endpoint's wall's lower boundaries
    if isinstance(pCell.below,float):
        pyMin=pCell.below
    elif isinstance(pCell.below,LineSegment):
        pyMin=pCell.below.m*p.x+pCell.below.b
    #check end endpoint's wall's lower boundaries
    if isinstance(qCell.below,float):
        qyMin=qCell.below
    elif isinstance(qCell.below,LineSegment):
        qyMin=qCell.below.m*q.x+qCell.below.b

    #check left endpoint's wall's upper boundaries
    if isinstance(pCell.above,float):
        pyMax=pCell.above
    elif isinstance(pCell.above,LineSegment):
        pyMax=pCell.above.m*p.x+pCell.above.b
    #check right endpoint's wall's upper boundaries
    if isinstance(qCell.above,float):
        qyMax=qCell.above
    elif isinstance(qCell.above,LineSegment):
        qyMax=qCell.above.m*q.x+qCell.above.b

    #create two walls
    return Wall(p.x,p.y,pyMin,pyMax),Wall(q.x,q.y,qyMin,qyMax)

"""
trimWalls: trim back the wall(s) if the new line segment has cutted through it
    Parameter(s):
        ls: the line segment
        walls:  a list of current walls
    Return:
        a list of walls that have been trimmed back
"""
def trimWalls(ls,walls):
    modifiedWalls=[]
    for i in range(len(walls)):
        currWall=walls[i]
        if currWall.x>ls.start.x and currWall.x<ls.end.x:
            currY=ls.m*currWall.x+ls.b
            if currY>currWall.ys and currY<currWall.ya:   #trim above
                currWall.ya=currY
                modifiedWalls.append(currWall)
            elif currY<currWall.ys and currY>currWall.yi: #trim below
                currWall.yi=currY
                modifiedWalls.append(currWall)
    return modifiedWalls
            
"""
buildTrapezoidalMap: uses randomized incremental algorithm to build the trapezoidal map
    Parameter(s)
        lineSegments:   a list of line segements
        dag:    the root of DAG
        cells:  a list of valid cells
"""
def buildTrapezoidalMap(lineSegments,dag,cells):
    walls=[]
    for ls in lineSegments:
        dag.printRecursive(0)
        p=ls.start
        q=ls.end
        #find the cells p and q are located
        pNode=pointLocate(p,dag)
        qNode=pointLocate(q,dag)
        pCell=pNode.itself
        qCell=qNode.itself
        #create and add walls
        pWall,qWall=wallCreate(pCell,qCell,p,q)
        walls.append(pWall)
        walls.append(qWall)
        walls=sorted(walls) #sort the walls based on their x values

        #case #2: both endpoints are inside the same cell
        if pCell==qCell:
            #remove the current cell
            cells.remove(pCell)
            #create new cells
            aCell=Cell(pCell.left,pWall,pCell.above,pCell.below)
            bCell=Cell(pWall,qWall,pCell.above,ls)
            cCell=Cell(pWall,qWall,ls,pCell.below)
            dCell=Cell(qWall,pCell.right,pCell.above,pCell.below)
            cells.append(aCell)
            cells.append(bCell)
            cells.append(cCell)
            cells.append(dCell)
            
            #update maps
            sNewNode=DAG(DAG(None,bCell,None),ls,DAG(None,cCell,None))
            pNode.itself=p
            pNode.leftChild=DAG(None,aCell,None)
            pNode.rightChild=DAG(sNewNode,q,DAG(None,dCell,None))
        
        #case 1&3: the line segment sets across different cells
        else:
            # Determine cells that are intersected by the line
            intersectedCells = [pCell]

            nextCell = pCell
            while nextCell != qCell:
                nextCandidates = [c for c in cells if (isinstance(c.left, Wall) and c.left.x == nextCell.right.x) or \
                                  (isinstance(c.left, float) and c.left == nextCell.right.x)]
                for c in nextCandidates:
                    if c.above.get_y(nextCell.right.x) > ls.get_y(nextCell.right.x) > c.below.get_y(nextCell.right.x):
                        nextCell = c
                        intersectedCells.append(c)
                        break

            prevTop = None
            prevBottom = None

            for c in intersectedCells:
                if c == pCell:
                    # Case 1: Left endpoint
                    leftCell = Cell(pCell.left, pWall, pCell.above, pCell.below)
                    topCell = Cell(pWall, pCell.right, pCell.above, ls)
                    bottomCell = Cell(pWall, pCell.right, ls, pCell.below)

                    cells.remove(pCell)
                    cells.append(leftCell)
                    cells.append(topCell)
                    cells.append(bottomCell)

                    pNode.itself = p
                    pNode.leftChild = DAG(None, leftCell, None)
                    pNode.rightChild = DAG(DAG(None, topCell, None), ls, DAG(None, bottomCell, None))

                    prevTop = pNode.rightChild.leftChild
                    prevBottom = pNode.rightChild.rightChild
                elif c == qCell:
                    # Case 1: Right endpoint
                    leftWall = c.left
                    if ls.get_y(leftWall.x) < leftWall.ys:
                        # Point is above new line, need to merge cell below

                        # Update the wall
                        qCell.left.yi = ls.get_y(qCell.left.x)

                        topCell = Cell(qCell.left, qWall, qCell.above, ls)
                        bottomCell = Cell(prevBottom.itself.left, qWall, ls, qCell.below)
                        rightCell = Cell(qWall, qCell.right, qCell.above, qCell.below)

                        cells.remove(qCell)
                        cells.append(rightCell)
                        cells.append(topCell)
                        cells.append(bottomCell)

                        qNode.itself = q
                        qNode.rightChild = DAG(None, rightCell, None)
                        qNode.leftChild = DAG(DAG(None, topCell, None), ls, DAG(None, bottomCell, None))

                        # Update the previous bottom node to be the same as this bottom node
                        cells.remove(prevBottom.itself)
                        prevBottom.itself = bottomCell
                    else:
                        # Point is below new line, need to merge cell above

                        # Update the wall
                        qCell.left.ya = ls.get_y(qCell.left.x)

                        topCell = Cell(prevTop.itself.left, qWall, qCell.above, ls)
                        bottomCell = Cell(qCell.left, qWall, ls, qCell.below)
                        rightCell = Cell(qWall, qCell.right, qCell.above, qCell.below)

                        cells.remove(qCell)
                        cells.append(rightCell)
                        cells.append(topCell)
                        cells.append(bottomCell)

                        qNode.itself = q
                        qNode.rightChild = DAG(None, rightCell, None)
                        qNode.leftChild = DAG(DAG(None, topCell, None), ls, DAG(None, bottomCell, None))

                        # Update the previous top node to be the same as this top node
                        cells.remove(prevTop.itself)
                        prevTop.itself = topCell
                else:
                    # Case 3: Cuts through trapezoid
                    leftWall = c.left
                    if ls.get_y(leftWall.x) < leftWall.ys:
                        # Point is above new line, need to merge cell below

                        # Update the wall
                        c.left.yi = ls.get_y(c.left.x)

                        # Create the two new cells
                        topCell = Cell(c.left, c.right, c.above, ls)
                        bottomCell = Cell(prevBottom.itself.left, c.right, ls, c.below)

                        cells.remove(c)
                        cells.append(topCell)
                        cells.append(bottomCell)


                        # Test a point to find this cell in the DAG
                        testX = (c.left.x + c.right.x) / 2
                        testPoint = Point(testX, ls.get_y(testX), None)
                        node = pointLocate(testPoint, dag)

                        # Update the node
                        node.itself = ls
                        node.leftChild = DAG(None, topCell, None)
                        node.rightChild = DAG(None, bottomCell, None)

                        # Update the previous bottom node to be the same as this bottom node
                        cells.remove(prevBottom.itself)
                        prevBottom.itself = bottomCell


                        prevTop = node.leftChild
                        prevBottom = node.rightChild
                    else:
                        # Point is below new line, need to merge cell above

                        # Update the wall
                        c.left.ya = ls.get_y(c.left.x)

                        # Create the two new cells
                        topCell = Cell(prevTop.itself.left, c.right, c.above, ls)
                        bottomCell = Cell(c.left, c.right, ls, c.below)

                        cells.remove(c)
                        cells.append(topCell)
                        cells.append(bottomCell)

                        # Test a point to find this cell in the DAG
                        testX = (c.left.x + c.right.x) / 2
                        testPoint = Point(testX, ls.get_y(testX), None)
                        node = pointLocate(testPoint, dag)

                        # Update the node
                        node.itself = ls
                        node.leftChild = DAG(None, topCell, None)
                        node.rightChild = DAG(None, bottomCell, None)

                        # Update the previous top node to be the same as this top node
                        cells.remove(prevTop.itself)
                        prevTop.itself = topCell

                        prevTop = node.leftChild
                        prevBottom = node.rightChild


def main():
    inputfile=open("D:/2020fall/716 computational geometry/coding assignments/a3/TrapezoidalMap/4.txt",'r')
    inputfile.readline()    #ignore number of line segments

    #get bounding box
    line=inputfile.readline().split()
    xMin=float(line[0])
    yMin=float(line[1])
    xMax=float(line[2])
    yMax=float(line[3])

    #intialize trapezoidal map
    cells=[]
    #initialCell=Cell(xMin,xMax,yMax,yMin)
    initialCell = Cell(Wall(xMin, yMax, yMin, yMax), Wall(xMax, yMin, yMin, yMax), LineSegment(Point(xMin, yMax, 'P'), Point(xMax, yMax, 'Q')), LineSegment(Point(xMin, yMin, 'P'), Point(xMax, yMin, 'Q')))
    cells.append(initialCell)
    dag=DAG(None,initialCell,None)

    #get points and line segments
    points=set()
    lineSegments=[]
    while True:
        line=inputfile.readline().split()
        if not line:
            break
        if float(line[0])<=float(line[2]):
            start=Point(float(line[0]),float(line[1]),'P')
            end=Point(float(line[2]),float(line[3]),'Q')
        else:
            start=Point(float(line[2]),float(line[3]),'P')
            end=Point(float(line[0]),float(line[1]),'Q')
        points.add(start)
        points.add(end)
        lineSegments.append(LineSegment(start,end))
    inputfile.close()

    plt.axis([xMin, xMax, yMin, yMax])

    buildTrapezoidalMap(lineSegments,dag,cells)

    dag.printRecursive(0)
    dag.plot("D:/2020fall/716 computational geometry/coding assignments/a3/TrapezoidalMap/plot.png")
    #TODO print out the required matrix

    #TODO ask user to input a point and call pointLocate to find its cell

if __name__ == "__main__":
    main()