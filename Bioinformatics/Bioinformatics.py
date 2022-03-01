from typing import List
from Blosumfile import blosum
import copy


List=['P','N','T','C','G','G','G','G','G','G','M','P','E','R','P']

TrainingList=['P','N','T','C','G','G','G','G','G','G','M','P','E','R','P']
#TrainingList="PNTCGGGGGGMPERP"
##########################################################Step 1 #############################################################################
for counter in range(0,len(List)-1):    
    
    if List[counter]==List[counter+1]:
        TrainingList[counter]='X'
        TrainingList[counter+1]='X'

 ##########################################################Step 2 #############################################################################
WordList=[]
WordLength=3
for Counter in range(0,len(List)-3):
    xList=[TrainingList[Counter],TrainingList[Counter+1],TrainingList[Counter+2]]
    WordList.append(xList)
       
##########################################################Step 3 #############################################################################

AminoAcids=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']

TheOutPutList = []

for word_index,word in enumerate(WordList):
    for litter_index,litter in enumerate(word):
        for Acid in AminoAcids :
            newlist = copy.deepcopy(word)
            newlist[litter_index] = Acid
            TheOutPutList.append(newlist)
    
           
           

##########################################################Step 4 #############################################################################
seeds=[]
Counter=0  
WhileCounter=60

for item in range(len(WordList)):
    
    while Counter <WhileCounter:
     value=0
     for litter in range(0,3): 
       if WordList[item][litter]==TheOutPutList[Counter][litter]:
          continue 
       try:
        value+=blosum[(WordList[item][litter],TheOutPutList[Counter][litter])]
       except:
             continue
     if value>13: 
      seeds.append(TheOutPutList[Counter])
     Counter+=1
    WhileCounter+=60


for OuterCounter in range(0,len(WordList)-1):
    value=0
    for InnerCounter in range(0,WordLength):
        try:
         value+=blosum[(WordList[OuterCounter][InnerCounter],WordList[OuterCounter][InnerCounter])]
        except:
             continue
    if value>13: 
         seeds.append(WordList[OuterCounter])

##########################################################step 5#########################################################################################################################################################

DatabaseSequence="RDQAHIFSLPEFFVPNTMPTQQMNKWKHILMSDDDDSARRCKGFYYWSFHGVNQQQGGLQEGGAQNTCHFIHHKGLPSTWSDQIYVFLHGYYCREGVPSASSPQRLHPQPISDICKSECHEHEYGAYPDRPRAMIFMFMFPPWKHILMPEMLPESDDDDSARQAHIFSFFVMPTQGVPSAS"
pointer=[]
QueryPointer=[]
for seedscounter in range(len(seeds)):#3
    for counter in range(len(DatabaseSequence)-3):  #178
      if(DatabaseSequence[counter]==seeds[seedscounter][0] and DatabaseSequence[counter+1]==seeds[seedscounter][1] and DatabaseSequence[counter+2]==seeds[seedscounter][2]):   
          pointer.append(counter)
          pointer.append(counter+2)
         
for seedscounter in range(len(seeds)):#3
    for counter in range(len(TrainingList)-3):  #178
      if(TrainingList[counter]==seeds[seedscounter][0] and TrainingList[counter+1]==seeds[seedscounter][1] and TrainingList[counter+2]==seeds[seedscounter][2]):   
          QueryPointer.append(counter)
          QueryPointer.append(counter+2)

frontQueryPointer=[]
BackQueryPointer=[]
frontDatabasePointer=[]
BackDatabasePointer=[]
for counter in range(0,len(QueryPointer)):
  if(counter%2==0):
   frontQueryPointer.append(QueryPointer[counter])
  else:
   BackQueryPointer.append(QueryPointer[counter])

for counter in range(0,len(pointer)):
  if(counter%2==0):
   frontDatabasePointer.append(pointer[counter])
  else:
   BackDatabasePointer.append(pointer[counter])

##################################################################step 6###########################################################################
def get_litter_score(db_index,trian_index):
    
    db_value = DatabaseSequence[db_index] # 14 P
    train_value = TrainingList[trian_index] # 0 P
    try :
        value = blosum[(train_value,db_value)]
        return value
    except :
        return 0  
    
def calculate_sequence_score(db_front,db_back,train_front,train_back):#14 -17 0-3

    full_score = 0                            # 14 15 16 17                             0 1 2 3
    for db_litter_index,train_litter_index in zip(range(db_front,db_back+1),range(train_front,train_back+1)):
        
        score = get_litter_score(db_litter_index,train_litter_index)
        # print("litter db index is {} , train index is {} , score is {}".format(db_litter_index,train_litter_index,score))
        full_score = full_score + score
    
    return full_score

def define_sequnce_and_the_extented_letters(db_front,db_back,train_front,train_back,startLenght,stard_db_index): # 14 17 0 3 #3 14
    #input   14 17 0 3 - 3 14 output "PNTQ"
    main_sequence = []
    for i in range(stard_db_index,stard_db_index+startLenght ): # 14 15 16
        main_sequence.append(i) # [14,15,16]

    seq = ""
    for i,tr in zip(range(db_front,db_back+1 ),range(train_front,train_back+1)):
        if i not in main_sequence : # 10 not in [14,15,16] PNT
            seq = seq + DatabaseSequence[i] + TrainingList[tr]
        else :
            seq = seq + DatabaseSequence[i]



    return seq

def if_extenction_is_good(db_front,db_back,train_front,train_back):
    if 0 <= db_front -1  < len(DatabaseSequence) and 0 <=  train_front -1 < len(TrainingList) :
        db_front = db_front - 1
        train_front = train_front - 1
    
    
    if 0 <= db_back +1  < len(DatabaseSequence) and 0 <=  train_back +1 < len(TrainingList) :
        db_back = db_back + 1
        train_back = train_back + 1
    
    return (db_front,db_back,train_front,train_back)

def extined_the_seqence(db_front,db_back,train_front,train_back,t,x): # 14 16 0 2 
    startLenght = (db_back - db_front)+1 # 3
    stard_db_index = db_front # 14 

    db_front,db_back,train_front,train_back = if_extenction_is_good(db_front,db_back,train_front,train_back)# 14 17 0 3

    topscore = calculate_sequence_score(db_front,db_back,train_front,train_back)# 14 17 0 3 

    if t > topscore :
        return None 
                # last_seq = (PNTQ , 15)
    last_seq = (define_sequnce_and_the_extented_letters(db_front,db_back,train_front,train_back,startLenght,stard_db_index) , topscore)
    
    while(True):

        db_front,db_back,train_front,train_back = if_extenction_is_good(db_front,db_back,train_front,train_back)

        score = calculate_sequence_score(db_front,db_back,train_front,train_back)

        if t > score :
            return last_seq

        if score > topscore :
            continue
        else :
            res =  topscore - score
           
            if res > x :
               return last_seq
            else:
                last_seq = (define_sequnce_and_the_extented_letters(db_front,db_back,train_front,train_back,startLenght,stard_db_index) , score)



#####################################################################step 7##########################################################################     
def extined_all_the_hits(frontDatabasePointer,BackDatabasePointer,frontQueryPointer,BackQueryPointer,t,x):
    Sequence_map  = []
    for db_front,db_back,train_front,train_back  in zip(frontDatabasePointer,BackDatabasePointer,frontQueryPointer,BackQueryPointer):
        res =  extined_the_seqence(db_front,db_back,train_front,train_back,t,x)
        Sequence_map.append(res)
    
    return Sequence_map


print(extined_all_the_hits(frontDatabasePointer,BackDatabasePointer,frontQueryPointer,BackQueryPointer,13,3))

        


# extined_all(frontDatabasePointer,BackDatabasePointer,frontQueryPointer,BackQueryPointer,t,x)
# [14, 65, 147]
# [16, 67, 149]
# [0, 1, 10]
# [2, 3, 12]
# 13
# 2

# (14,16,0,2) 13 2 --> extined_one_seqence()
# (65, 67, 1, 3) 13 2 -->extined_one_seqence()
# (147, 149, 10, 12) 13 2 --> extined_one_seqence()


# 14,16,0,2

# database 14  17
# TrainingList  3

# 14 17 0 3

# PNT --> PNTQ