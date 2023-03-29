from calendar import c
import random
import copy
import math
import matplotlib.pyplot as plt

class Structure:
    def __init__(self, seqn):
        self.seq = seqn
        self.buriedScore = 0
        self.pairwiseScore = 0
        self.sizeOfMatrix = (len(seqn) * 2) + 1
        self.matrix = [[[9 for x in range(self.sizeOfMatrix)]
                       for y in range(self.sizeOfMatrix)] for z in range(self.sizeOfMatrix)]
        self.moveSet = []
        self.score = -1
        # Direction Map ()==> Up=0, Down=1, Left=2, Right=3, 3D-Up=4, 3D-Down=5)
        # Move Change in 3D Lattice
        self.dir = {0: (-1, 0, 0), 1: (1, 0, 0), 2: (0, -1, 0), 3: (0, 1, 0), 4: (0, 0, 1), 5: (0, 0, -1)}
        self.pointList = [0]*len(seqn)
        self.bondingPoints = [[[],[],[]],[[],[],[]],[[],[],[]],[[],[],[]]]
        self.__initMoveSet()

    # First Random Structure and MoveSet
    def __initMoveSet(self):
        i, j, k = len(self.seq), len(self.seq), len(self.seq)
        self.matrix[i][j][k] = self.seq[0]
        self.pointList[0] = (i,j,k)
        c = 1
        z=1
        while c < len(self.seq):
            # Chosing random next move
            curr = random.randint(0, 5)
            # Directions Up=0, Down=1, Left=2, Right=3
            ti = i + self.dir[curr][0]
            tj = j + self.dir[curr][1]
            tk = k + self.dir[curr][2]
            if self.matrix[ti][tj][tk] == 9:
                self.matrix[ti][tj][tk] = self.seq[c]
                self.moveSet.append(curr)
                i = ti
                j = tj
                k = tk
                self.pointList[z] = (i,j,k)
                z+=1
                c += 1
            # If crosswalk occurs
            else:
                ti = i
                tj = j
                tk = k

    # New Structure after Mutation
    def resetMatrix(self, mSet):
        self.matrix = [[[9 for x in range(self.sizeOfMatrix)]
                       for y in range(self.sizeOfMatrix)] for z in range(self.sizeOfMatrix)]
        # Starting Index of sequence
        i, j, k = len(self.seq), len(self.seq), len(self.seq)
        self.matrix[i][j][k] = self.seq[0]
        self.pointList[0] = (i,j,k)
        flag = False
        z = 1
        # Drawing Structure against moveset
        for x in range(0, len(mSet)):
            mov = mSet[x]
            i += self.dir[mov][0]
            j += self.dir[mov][1]
            k += self.dir[mov][2]
            # If crosswalk happens, structure is not feasible
            if self.matrix[i][j][k] == 9:
                self.matrix[i][j][k] = self.seq[x+1]
                self.pointList[z] = (i,j,k)
                z += 1
            else:
                flag = True
                break
        return flag

    # Given Bond score and All Possible Bonds,Buried Bonds
    def getBondScore(self, i, j,k, p, q,r):
        # If Bond is P-P (Returns Buried or Pairwise)
        if self.matrix[i][j][k] == self.matrix[p][q][r] == 0:
            return -2, 0
        # If Bond is H-H
        elif self.matrix[i][j][k] == self.matrix[p][q][r] == 1:
            return -3, 0
        # If Bond is P-H
        elif self.matrix[i][j][k] == 0 and self.matrix[p][q][r] == 1:
            return 1, 0
        # If Bond is H-P
        elif self.matrix[i][j][k] == 1 and self.matrix[p][q][r] == 0:
            return 1, 0
        # If H has Buried Neighbour
        elif self.matrix[i][j][k] == 1 and self.matrix[p][q][r] == 9:
            return 0, 1
        else:
            return 0, 0

    # Neighbouring Bond, Burried (in case of H Score)
    def calNodeScore(self, i, j, k, ms):
        # curr moveSet action
        p, q, r = -2, -2, -2
        if ms < len(self.moveSet):
            p = i + self.dir[self.moveSet[ms]][0]
            q = j + self.dir[self.moveSet[ms]][1]
            r = k + self.dir[self.moveSet[ms]][2]
        # Prev moveSet action
        x, y, z = -2, -2, -2
        # If MoveSet index is not at 0, If at 0 we can't access 0-1
        if ms - 1 >= 0:
            move = self.moveSet[ms-1]
            # If Move is from UP to DOWN, Now Standing at DOWN previous is at UP
            # Same for RIGHT and LEFT
            if move == 0:
                move = 1
            elif move == 1:
                move = 0
            elif move == 2:
                move = 3
            elif move == 3:
                move = 2
            elif move ==4:
                move = 5
            else:
                move = 4
            x = i + self.dir[move][0]
            y = j + self.dir[move][1]
            z = k + self.dir[move][2]
        s = len(self.seq) * 2
        pairwise = 0
        buried = 0
        # Visiting Neighbour Indexes
        for n in range(i-1, i+2):
            for m in range(j-1, j+2):
                for l in range(k-1, k+2):
                # Checking if the Index is a Node in Alternative in Sequence
                    if n == i and m == j and l == k or n == p and m == q and l == r or n == x and m == y and l == z:
                        continue
                    else:
                        pw, br = self.getBondScore(i, j,k, n, m,l)
                        pairwise += pw
                        buried += br
        return pairwise, buried

    # Applying Energy Function
    def computeEnergy(self):
        # index of first move
        i, j, k = len(self.seq), len(self.seq), len(self.seq)
        # first residue score where its move is at index 0
        pairwise, buried = self.calNodeScore(i, j, k, 0)
        for x in range(0, len(self.moveSet)):
            # Move in MoveSet
            mov = self.moveSet[x]
            i += self.dir[mov][0]
            j += self.dir[mov][1]
            k += self.dir[mov][2]
            # Single(not first) H/P Neighbouring Score Counting
            pw, br = self.calNodeScore(i, j, k, x+1)
            pairwise += pw
            buried += br
        pairwise = pairwise // 2
        # Energy Function
        self.score = 2*buried + 3*pairwise
        self.buriedScore = buried
        self.pairwiseScore = pairwise
        return self.score

    # Converting Singular point List to Double Points as Line Start and End
    def pointToLine(self):
        tempPoint = [0]*(len(self.pointList)-1)
        for i in range(0, len(self.pointList)-1):
            tempPoint[i] = [self.pointList[i], self.pointList[i+1]]
        return tempPoint


    # Displays Best Structure along with score
    def printMatrix(self):
        fig = plt.figure(figsize=(10,10))
        x,y,z = [0]*len(self.pointList), [0]*len(self.pointList), [0]*len(self.pointList)
        ax = fig.add_subplot(111, projection='3d')
        for i in range(0,len(self.pointList)):
            x[i] = self.pointList[i][0]
            y[i] = self.pointList[i][1]
            z[i] = self.pointList[i][2]
            
        for i in range(0,len(z)):
            if i==0:
                ax.scatter(x[i],y[i],z[i], c= '#e74c3c',s=500.0)
            elif self.seq[i] == 1:
                ax.scatter(x[i],y[i],z[i], c= '#2ecc71',s=500.0)
            else: 
                ax.scatter(x[i],y[i],z[i], c= '#3498db',s=500.0)
        ax.plot(self.bondingPoints[0][0],self.bondingPoints[0][1],self.bondingPoints[0][2])
        ax.plot(x,y,z,c='#000000',linewidth=5)
        plt.show()

# Changes to new move set with certain restrictions to converge in less time
def mutation(struc: Structure, count):
    flag = True
    for a in range(count):
        while flag:
            mSet = copy.deepcopy(struc.moveSet)
            # Index of MoveSet where the change will happen
            i = random.randint(1, (len(struc.moveSet)-1))
            m = random.randint(0, 5)
            # New Move is equal to itself
            if mSet[i] == m:
                continue
            elif i == 0:
                mSet[i] = m
                flag = struc.resetMatrix(mSet)
            # Change should not move to itself
            elif mSet[i-1] == 0 and m == 1:
                continue
            elif mSet[i-1] == 1 and m == 0:
                continue
            elif mSet[i-1] == 2 and m == 3:
                continue
            elif mSet[i-1] == 3 and m == 2:
                continue
            elif mSet[i-1] == 4 and m == 5:
                continue
            elif mSet[i-1] == 5 and m == 4:
                continue
            else:
                mSet[i] = m
            # Structure got from new move set and structure has no crosswalk
                flag = struc.resetMatrix(mSet)
        struc.moveSet[i] = m

# Always Accepting Better moveSet but Also accepting worse with probability
def simAnnealing(count, hpSeq):
    seq = hpSeq
    s1 = Structure(seqn=seq)
    s3 = copy.deepcopy(s1)
    e3 = s3.computeEnergy()
    # Initial Temperature tends to Low from Current
    temp = count
    t = temp
    for i in range(count):
        s2 = copy.deepcopy(s1)
        mutation(s2, random.randint(1, 2))
        # Nest BEST/WORSE Enery
        e1 = s1.computeEnergy()
        e2 = s2.computeEnergy()
        # Accepting Better
        if e2 < e1:
            s1 = copy.deepcopy(s2)
        # Using Metropolis Acceptance Formula Where q will get Lower and chances of worse
        # solution selection as well
        else:
            v = (e2-e1)/(t)
            q = math.exp(-v)
            # q decreases grarually and less chances to accept worse solution
            if random.random() > q:
                continue
            else:
                s1 = copy.deepcopy(s2)
        # Keeping Best
        if e2 < e3:
            s3 = copy.deepcopy(s1)
            e3 = e2
        # T decreases
        t = temp - (i)
    # s1 is Solution of Last Iteration, s3 is Best Solution
    return s1, s3,

# Mapping H to 1 and P to 0
def sequenceMap(hpSeq):
    mappedhpSeq = []
    rightSeq = True
    if len(hpSeq) >= 3:
        for i in range(len(hpSeq)):
            if hpSeq[i] == "H" or hpSeq[i] == "h":
                mappedhpSeq.append(1)
            elif hpSeq[i] == "P" or hpSeq[i] == "p":
                mappedhpSeq.append(0)
            else:
                rightSeq = False
    else:
        rightSeq = False
    return mappedhpSeq, rightSeq

def main():
    runningFlag = True
    while runningFlag == True:
        hpSeqAlpha = input("Enter Sequence and Press ENTER!   ")
        hpSeqBeta, rightSeq = sequenceMap(hpSeq=hpSeqAlpha)
        if rightSeq == True:
            noOfIterations = input(
                "Enter NUMBER of Iterations and Press ENTER!   ")
            fs, bs = simAnnealing(int(noOfIterations), hpSeqBeta)
            print("BEST STRUCTURE SCORE:", end=" ")
            print(bs.score)
            print("BEST STRUCTURE BURIED SCORE:", end=" ")
            print(bs.buriedScore)
            print("BEST STRUCTURE PAIRWISE SCORE:", end=" ")
            print(bs.pairwiseScore)
            bs.printMatrix()
            runningFlag = False
        else:
            print(
                "Entered Sequence is INCORRECT, ENTER AGAIN! (Sequence should have H and P only and atleast 3 residues)")
            continue

if __name__ == "__main__":
    main()