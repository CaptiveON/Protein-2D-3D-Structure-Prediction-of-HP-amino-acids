import random
import copy
import math
import pygame

class Structure:
    def __init__(self, seqn):
        self.seq = seqn
        self.sizeOfMatrix = (len(seqn) * 2) + 1
        self.matrix = [[9 for x in range(self.sizeOfMatrix)]
                       for y in range(self.sizeOfMatrix)]
        self.moveSet = []
        self.score = -1
        # Direction Map ==> Up=0, Down=1, Left=2, Right=3 and Move Change in 2D Lattice
        self.dir = {0: (-1, 0), 1: (1, 0), 2: (0, -1), 3: (0, 1)}
        self.pointList = [0]*len(seqn)
        self.__initMoveSet()

    # First Random Structure and MoveSet
    def __initMoveSet(self):
        i, j = len(self.seq), len(self.seq)
        self.matrix[i][j] = self.seq[0]
        z = 1
        # self.pointList[0] = (i*50,j*50)
        c = 1
        while c < len(self.seq):
            curr = random.randint(0, 3)
            # Directions Up=0, Down=1, Left=2, Right=3
            ti = i + self.dir[curr][0]
            tj = j + self.dir[curr][1]
            if self.matrix[ti][tj] == 9:
                self.matrix[ti][tj] = self.seq[c]
                # self.pointList[z] = (ti*50,tj*50)
                z += 1
                self.moveSet.append(curr)
                i = ti
                j = tj
                c += 1
            else:
                ti = i
                tj = j

    # New Structure after Mutation
    def resetMatrix(self, mSet):
        self.matrix = [[9 for x in range(self.sizeOfMatrix)]
                       for y in range(self.sizeOfMatrix)]
        i, j = len(self.seq), len(self.seq)
        self.matrix[i][j] = self.seq[0]
        z = 1
        # self.pointList[0] = (i*50,j*50)
        flag = False
        for x in range(0, len(mSet)):
            mov = mSet[x]
            i += self.dir[mov][0]
            j += self.dir[mov][1]
            if self.matrix[i][j] == 9:
                self.matrix[i][j] = self.seq[x+1]
                # self.pointList[z] = (i*50,j*50)
                z += 1
            else:
                flag = True
                break
        return flag

    # Given Bond score and All Possible Bonds,Buried Bonds
    def getBondScore(self, i, j, p, q):
        # If Bond is P-P
        if self.matrix[i][j] == self.matrix[p][q] == 0:
            return -2, 0
        # If Bond is H-H
        elif self.matrix[i][j] == self.matrix[p][q] == 1:
            return -3, 0
        # If Bond is P-H
        elif self.matrix[i][j] == 0 and self.matrix[p][q] == 1:
            return 1, 0
        # If Bond is H-P
        elif self.matrix[i][j] == 1 and self.matrix[p][q] == 0:
            return 1, 0
        # If H has Buried Neighbour
        elif self.matrix[i][j] == 1 and self.matrix[p][q] == 9:
            return 0, 1
        else:
            return 0, 0

    # Neighbouring Bond, Burried (in case of H Score)
    def calNodeScore(self, i, j, ms):
        # curr moveSet action
        p, q = -2, -2
        if ms < len(self.moveSet):
            p = i + self.dir[self.moveSet[ms]][0]
            q = j + self.dir[self.moveSet[ms]][1]
        # Prev moveSet action
        x, y = -2, -2
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
            else:
                move = 2
            x = i + self.dir[move][0]
            y = j + self.dir[move][1]
        s = len(self.seq) * 2
        pairwise = 0
        buried = 0
        # Visiting Neighbour Indexes
        for n in range(i-1, i+2):
            for m in range(j-1, j+2):
                # Checking if the Index is a Node in Alternative in Sequence
                if n == i and m == j or n == p and m == q or n == x and m == y:
                    continue
                elif s > n >= 0 and s > m >= 0:
                    pw, br = self.getBondScore(i, j, n, m)
                    pairwise += pw
                    buried += br
        return pairwise, buried

    # Applying Energy Function
    def computeEnergy(self):
        i, j = len(self.seq), len(self.seq)
        pairwise, buried = self.calNodeScore(i, j, 0)
        for x in range(0, len(self.moveSet)):
            # Move in MoveSet
            mov = self.moveSet[x]
            i += self.dir[mov][0]
            j += self.dir[mov][1]
            # Single H/P Neighbouring Score Counting
            pw, br = self.calNodeScore(i, j, x+1)
            pairwise += pw
            buried += br
        pairwise = pairwise // 2
        # Energy Function
        self.score = 2*buried + 3*pairwise
        print("Score:", end=" ")
        print(self.score)
        return self.score

    # Converting Singular point List to Double Points as Line Start and End
    def pointToLine(self):
        tempPoint = [0]*(len(self.pointList)-1)
        for i in range(0, len(self.pointList)-1):
            tempPoint[i] = [self.pointList[i], self.pointList[i+1]]
        return tempPoint

    # Displays Best Structure along with score
    def printMatrix(self):
        screensize_width = 1250
        screensize_height = 800
        # If Sequence is Large, Lines will be smaller as used percentage of Sequence
        scaleFactor = 50*(1-((len(self.seq)*0.2)/100))
        li = screensize_height/2
        lj = screensize_height/2
        self.pointList[0] = (li, lj)
        curMove = 0
        # Mapping Moveset to PyGame Display
        for i in range(1, len(self.moveSet)+1):
            curMove = self.moveSet[i-1]
            li = li+(scaleFactor*self.dir[curMove][0])
            lj = lj+(scaleFactor*self.dir[curMove][1])
            self.pointList[i] = (li, lj)
        tempList = self.pointToLine()
        pygame.init()
        screen = pygame.display.set_mode((screensize_width, screensize_height))
        pygame.display.set_caption("H-P Model Protein Structure Prediction")
        # Images for H and P
        imageH1 = pygame.image.load('h.png')
        imageP0 = pygame.image.load('p.png')
        color = (0, 125, 255)
        bgColor = (255, 255, 255)
        startColor = (0, 255, 0)
        done = False
        screen.fill(bgColor)
        # Connectiog Bonds
        for i in range(0, len(tempList)):
            if i == 0:
                pygame.draw.lines(screen, startColor, False,
                                  tempList[i], width=3)
            else:
                pygame.draw.lines(screen, color, False, tempList[i], width=3)
        # Drawing H and P
        for i in range(0, len(self.pointList)):
            if self.seq[i] == 1:
                screen.blit(
                    imageH1, (self.pointList[i][0]-16, self.pointList[i][1]-16))
            else:
                screen.blit(
                    imageP0, (self.pointList[i][0]-16, self.pointList[i][1]-16))
        # Simple Side menu with Best Score and Length of Sequence
        pygame.draw.lines(screen, startColor, False, [
                          (850, 50), (850, 750)], width=3)
        font = pygame.font.Font('lato.ttf', 20)
        text5 = font.render(
            'Green Bond is Start of Sequence! ', True, color, startColor)
        textRect = text5.get_rect()
        textRect.center = (1065, 300)
        screen.blit(text5, textRect)
        text = font.render('Best Structure Energy: ', True, color, bgColor)
        textRect = text.get_rect()
        textRect.center = (1050, 400)
        screen.blit(text, textRect)
        text2 = font.render(str(self.score), True, startColor, bgColor)
        textRect = text2.get_rect()
        textRect.center = (1185, 400)
        screen.blit(text2, textRect)
        text3 = font.render('Length of Sequence: ', True, color, bgColor)
        textRect = text3.get_rect()
        textRect.center = (1050, 450)
        screen.blit(text3, textRect)
        text4 = font.render(str(len(self.seq)), True, startColor, bgColor)
        textRect = text4.get_rect()
        textRect.center = (1185, 450)
        screen.blit(text4, textRect)
        while not done:
            for event in pygame.event.get():
                if event.type == pygame.QUIT:
                    done = True
            pygame.display.flip()

# Changes to new move set with certain restrictions to converge in less time
def mutation(struc: Structure, count):
    flag = True
    for a in range(count):
        while flag:
            mSet = copy.deepcopy(struc.moveSet)
            # Index of MoveSet where the change will happen
            i = random.randint(1, (len(struc.moveSet)-1))
            m = random.randint(0, 3)
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
            elif [i-1] == 2 and m == 3:
                continue
            elif mSet[i-1] == 3 and m == 2:
                continue
            else:
                mSet[i] = m
            # Structure got from new move set
                flag = struc.resetMatrix(mSet)
        struc.moveSet[i] = m

# Mapping H to 1 and P to 0
def sequenceMap(hpSeq):
    mappedhpSeq = []
    rightSeq = True
    for i in range(len(hpSeq)):
        if hpSeq[i] == "H" or hpSeq[i] == "h":
            mappedhpSeq.append(1)
        elif hpSeq[i] == "P" or hpSeq[i] == "p":
            mappedhpSeq.append(0)
        else:
            rightSeq = False
    return mappedhpSeq, rightSeq

# Always Accepting Better moveSet but Also accepting worse with probability
def simAnnealing(count, hpSeq):
    seq = hpSeq
    s1 = Structure(seqn=seq)
    s3 = copy.deepcopy(s1)
    e3 = s3.computeEnergy()
    # Initial Temperature tends to Low from Current
    temp = 10
    t = temp
    for i in range(count):
        s2 = copy.deepcopy(s1)
        mutation(s2, random.randint(1, 2))
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
    t = temp - (i*1.2)
    # s1 is Solution of Last Iteration, s3 is Best Solution
    return s1, s3,

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
            bs.printMatrix()
            runningFlag = False
        else:
            print(
                "Entered Sequence is INCORRECT, ENTER AGAIN! (Sequence should have H and P only)")
            continue

if __name__ == "__main__":
    main()