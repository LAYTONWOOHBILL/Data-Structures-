class Point:

    """Summary of class here.

    Create a class Point with xcoord and ycoord

    Attributes:
        __repr__: A boolean indicating if we like SPAM or not.
        __gt__: Compare two point xcoord and ycoord, first on ycoord then xcoord
        __lt__: Compare two point xcoord and ycoord, first on ycoord then xcoord
        __eq__: Compare two point have same xcoord and ycoord
    """
    def __init__(self,xcoord=0,ycoord=0):
        # attributes xcoord and ycoord for a Point
        self.x = xcoord
        self.y = ycoord

    def __repr__(self):
        # representation of Point class objects
        return "Point({}, {})".format(self.x,self.y)


#readfile and add point into list
def addpoints_grahamAlgorithm(points,file):

    """readfile and add point into list

    Args:
        points: Empty list
        file: external file

    Returns:
        a pointlist with point been add
        example:
            [Point(0, 5), Point(-1, 3), Point(2, 4), Point(0, 0), Point(-1, -1), Point(2, 0), Point(1, -2)]

    Raises:
      if readline error

    """
    f = open(file, "r")                                 # open file
    line = f.readline()                                 # readline

    while line:                                         # while can read line do

        try:                                            # try
            line = line.split(" ")                          # ['0', '5\n']
            line[1] = line[1].rstrip("\n")                  # drop '\n'
            point = Point(eval(line[0]),eval(line[1]))      # cast str to int or float
            points.append(point)                            # append in list
            line = f.readline()                             # read next line
        except:                                         # raises error
            print("Error")                                  # print

    return


# Find a point with minimun y-coordinate then (minimun x-coordinate.)
def findmin(points):

    """Find the min point in Pointslist with minimun y-coordinate then ( minimun x-coordinate.)

    Args:
        points: Point list

    Returns:
        a min point in list
        example:
            Point(1, -2)

    """
    if len(points) == 0:                                # Empty lists
        return print ("no point in points list")

    else:
        minp = points[0]                                # set first one to min

        for i in points[1:]:                            # looping                                       EX: Min = Point(17, 2)  Points: Point(4, 1), Point(3, 1)
            if i.y < minp.y:                                # if that-ycoord less then min-ycoord       EX: Min = Point(17, 2) < Points: Point(4, 1)
                minp = i                                    # set new min point to that point           EX: Min = Point(4, 1)
            elif i.y == minp.y:                             # if that-ycoord eqaul then min-ycoord      EX: Min = Point(4, 1) = Points: Point(3, 1)
                if i.x < minp.x:                                # compare two x point if that-xcoord less then min-xcoord
                    minp = i                                    # set new min point to that point       EX: Min = Point(4, 1) > Points: Point(3, 1)
            else:                                           # if that-ycoord greater than min-ycoord    EX: Min = Point(3, 1)
               pass                                         # pass

        return minp


 # 2. Sort all the other point by comparing slope with the starting point

    """
        Original I try just use cross product, then polar angel
        but both doesnt help me sort all the points from right to left.
    """

def sortslope(points,starting_point):
    
    """Find all the points slope with the starting point

    slope = (point_i.ycoordinate-point_j.ycoordinate)/(point_i.xcoordinate-point_j.xcoordinate)

    Args:
        points: Point list
        starting_point: min Point

    Returns:
        a sorted point list by right to left from the starting_point

    """

    points.remove(starting_point)                                               # not finding starting_point slope by itself
    slope_list=[]                                                               # create an empty list for saving slope for sort later

    for point in points:                                                        # find each point slope with starting_point
        if (starting_point.x-point.x) != 0:                                         # if (starting_point.x-point.x) is not zero
            slope = (starting_point.y-point.y)/(starting_point.x-point.x)               # we can do divison noramlly
            slope_list.append((point,slope))                                            # append in list EX. [(Point(5, 4), 1.5), (Point(4, 3), 2.0)]
        else:                                                                       # else (starting_point.x-point.x) is zero
            slope_list.append((point,float('inf')))                                     # starting_point and this.point are on the same x-coordinate
    slope_list.sort(key=lambda slope: slope[1])                                 # sorted by slope  [(Point(6, 2), 0.5), (Point(5, 4), 3.0)]

    sortslopelist=[]                                                            # create an empty list for saving slope for sort later
    sortslopelist.append(starting_point)
                                                             
    for i in slope_list:                                                        # put positive slope in new list from min to max
        if i[1] >= 0:                                                           # EX.  (Point(2, 0), 2.0), (Point(2, 4), 6.0)
            sortslopelist.append(i[0])                                          # [Point(2, 0), Point(2, 4)] 
  
    for i in slope_list:                                                        # put negative slope in new list from min to max
        if i[1] < 0:                                                            # EX. (Point(0, 5), -7.0), (Point(-1, 3), -2.5)
            sortslopelist.append(i[0])                                          # [Point(2, 0), Point(2, 4), Point(0, 5), Point(-1, 3)]   

    points = sortslopelist
    
    return points

# vector
def vector(point_i,point_j):
     
    """compute a vector of two points

      vector = (point_j.x - point_i.x , point_j.y - point_i.y)

    Args:
        point_i = Point.this
        point_j = Point.other
        
    Returns:
        a vector of two points
    """
    
    return ((point_j.x - point_i.x , point_j.y - point_i.y))                   

# cross product
def cross_product (vector_i,vector_j):
    
    """compute a cross_product of two vectors

    cross_product = (vector_i.x * vector_j.y)-(vector_i.x *vector_j.y)

    Args:
        vector_i = vector_i.this
        vector_j = vector_j.other
        
    Returns:
         a cross_product of two vectors
    """
    return ((vector_i[0] * vector_j[1])-(vector_i[1] *vector_j[0]))

def compute_cross_product(point_top,point_next_top,point_i):
    
    """compute a cross_product of three points

    Args:
        point_top = stack top
        point_next_top = stack second top
        point_i = new point
        
    Returns:
         a cross_product of two three points
    """
    
    t_nt_vector = vector(point_top,point_next_top)                                  #return stack top and stack second top vector
    t_j = vector(point_top,point_i)                                                 #return stack top and new point vector
    return cross_product(t_nt_vector,t_j)                                           #return cross product


# initialise a stack Γ = ∅
def init_stack(points):
    
    init_stack = []                                                                 # create a empty stack

    if len(points) < 3:                                                             # if input point less than three 
        for point in points:                                                        # we just push points no need to compare with cross product
            init_stack.append(point)
            return init_stack
    else:                                                                           # else
        init_stack.append(points[0])                                                    # push(Γ, p1); push(Γ, p2);
        init_stack.append(points[1])

        i = 2;                                                                          # init index for later starting.
        while len(init_stack) < 3 and i < len(points):                                  # while stack has not contain 3 point yet or we are not run out of the point in points
            point_i= points[i]                                                              # set point_i = new point which we are trying to put in 
            top = init_stack[-1]                                                            # stack top
            next_top = init_stack[-2]                                                       # stack second top
            cp_value= compute_cross_product(top,next_top,point_i)                           # compute three points cross product
            if cp_value == 0:                                                               # if three points are collinear
                init_stack.pop()                                                                # we pop the top one
            init_stack.append(points[i])                                                    # put the new points on it since its collinear
            i+=1                                                                        
        return init_stack, points.index(init_stack[-1])+1                            # return stack and index for starting convex hull begain point


def GrahamAlgorithm(file):
    
    """Graham’s Algorithm: 

    1. Find a point in S having minimun y-coordinate with minimun x-coordinate. 

    2. Sort all the other point by Polar angle  // Operation: Comparison Runtime: O(1)
    // Compare P1 to Pi and P1 to Pj, we can use cross product;
    // Pi < Pj iff <P1,Pi, Pj> is a left turn;
    let the sorted list be ⟨p2, . . . , pn⟩      

    3. initialise a stack Γ = ∅;

    4. push(Γ, p1); push(Γ, p2); push(Γ, p3);

    5.
        for i = 4 to n do;   //Runtime: O(n)
            while ⟨NextToTop(Γ),Top(Γ), pi⟩ is not a left turn do //Runtime: O(less than n)
            PoP(Γ);
    Push(Γ, pi);

        // Once a point is popped from Γ it will never be reinserted
        // never pop same point more than once --> it will never take more than n
        // O(n) * O(less n) = O(n)*C
    6.
        Output the points in Γ as the vertices of CH(S);

    """

    points = list()                                                             # create an empty list for puting point
    addpoints_grahamAlgorithm(points,file)                                      # add points by reading txt
    starting_point = findmin(points)                                            # Find a point in S having minimun y-coordinate with minimun x-coordinate
    points=sortslope(points,starting_point)                                     # sort point by slope
    stack , starting_index = init_stack(points)                                 # initialise a stack and push three point in stack first
    
    for point_i in points[starting_index:]:                                     # for each point in Points where start from begain point after we do the initialise a stack
        top = stack[-1]                                                         # stack top
        next_top = stack[-2]                                                    # stack second top
        cp_value = compute_cross_product(top,next_top,point_i)                  # compute three points cross product
        while cp_value >= 0 and len(stack)>0:                                       # while ccross product right turn and stack is not empty do
            stack.pop()                                                             # we pop the top one
            top = stack[-1]                                                         # new stack top
            next_top = stack[-2]                                                    # stack second top
            cp_value = compute_cross_product(top,next_top,point_i)                  # compute three points cross product again
        stack.append(point_i)                                                   # left turn, push in the stack 

    return stack


def iteration_string(stack):
    
    """print out nice format for points

    Args:
        stack = convex_hull
        
    """
    
    it_str = ""

    for point in stack:
        it_str=it_str+str(point)+' '
    print(it_str)


def main():
    
    print("#Author: Wilson Wu\n#Date:2021.04.12\n#StudentID:1939700\n")
    
    file_list = ["case1.txt","case2.txt","case3.txt","case4.txt","case5.txt","case6.txt"]
    
    for file in file_list:
        print('='*120)
        print("Filename: "+file)
        convex_hull=GrahamAlgorithm(file)
        iteration_string(convex_hull)
        print('='*120+'\n')
        
        
if __name__ == "__main__":
    main()
