#!/usr/bin/env python2
# -*- coding: utf-8 -*-
__prj__ = 'Motivational Go/NoGo Task'
__license__ = 'See attached licence'
__author__ = 'Johannes Algermissen'
__email__ = 'j.algermissen@donders.ru.nl'
__date__ = '2019/10/26'

###########################################################################
# Load modules
###########################################################################

# my own modules
import os.path # for absolute path
from psychopy import core, visual, event, data, gui
from Tkinter import *
# import pygame
import random

# track logging
# import logging
# logger = logging.getLogger()
# logger.setLevel(logging.DEBUG)
###########################################################################
# Request participant number, gender, and age via GUI
###########################################################################

def function_gui_info():
    guiDict = {"A. Subject":999,"B. Age":999,"C. Gender":["male","female","other"],"D. Hand":["leftHand","rightHand"],
    "E. Practice":["1: Yes","2: No"], "F. Block 1":["1: Yes","2: No"], "G. Block 2":["1: Yes","2: No"], "H. Block 3":["1: Yes","2: No"], "I. Block 4":["1: Yes","2: No"], "J. nRep":20} # define dictionary with infos to be asked for
    dialog = gui.DlgFromDict(dictionary = guiDict, title = "Go/NoGo Task", fixed = ("options")) # ask for infos
    if dialog.OK: # if info received: define variables
        subject = guiDict["A. Subject"]
        age = int(guiDict["B. Age"])
        gender = guiDict["C. Gender"]
        hand = guiDict["D. Hand"]
        if guiDict["E. Practice"] == "1: Yes":
            isPractice = True
        else: 
            isPractice = False
        if guiDict["F. Block 1"] == "1: Yes":
            isBlock1 = True
        else: 
            isBlock1 = False
        if guiDict["G. Block 2"] == "1: Yes":
            isBlock2 = True
        else: 
            isBlock2 = False
        if guiDict["H. Block 3"] == "1: Yes":
            isBlock3 = True
        else: 
            isBlock3 = False
        if guiDict["I. Block 4"] == "1: Yes":
            isBlock4 = True
        else: 
            isBlock4 = False
        nTrials = int(guiDict["J. nRep"]) * 16; # 4 blocks are 4 stimuli, so 16*nRep stimuli in total
        exists = os.path.isfile("TrialHandler/stimuluslist_test_stimNum_{}_sub_{:0>3d}_Part1.csv".format(nTrials,subject))        
        if not exists:
            print(u"Stop -- trialHandler file for subject {} with {} trials does not exist".format(subject,nTrials))
            core.quit()
        return (subject, age, gender, hand, isPractice, isBlock1, isBlock2, isBlock3, isBlock4, nTrials)
    else: # otherwise quit experiment and display reason
        print("Stop -- missing participant information")
        core.quit()
        

###########################################################################
# Waitkeys in pygame
###########################################################################

def waitforbutton(waitTime):
    core.wait(waitTime) # make sure that participants see instructions for a certain amount of time, and do not just click through it        
    keyPress = event.waitKeys(keyList = ['space','escape'], maxWait = 600) # maximally wait 600 seconds, i.e. 10 minutes (in case subject gets experimenter to show them something)
#    print(u'1. Pressed key is: {}'.format(str(keyPress[0]))) # Print pressed key
    if keyPress[0] == 'escape': # check if keypress was abort
        core.quit()

###########################################################################
# Overall function displaying slide of instructions
###########################################################################

def function_display_instructions(file, waitTime = 0.1): # default 1.5 sec
    instr_image.setImage(file) # retrieve new image
    instr_image.draw()
    win.flip()
    waitforbutton(waitTime)
    win.flip() # clean screen 
    core.wait(0.3) # wait for 300 ms until next screen follows (avoid rapid transitions)
    
###########################################################################
# Display countdown at beginning of each block
###########################################################################

def function_countdown():
    for i in ["3","2","1","Start!"]:
        if 'escape' in event.getKeys(): # exit with Escape
            core.quit()
        center_fix.setText(i)
        center_fix.draw()
        win.flip() 
        core.wait(1) # make count down slower
    
###########################################################################
# Run trials -- either for practice or for test
###########################################################################

def function_run_trials(trials, block, part):
    # Display countdown
    function_countdown()
    # Initial fixation cross to calm down pupil:
    center_fix.setText('+')
    center_fix.draw()
    win.flip()
    # Intialize points to-be-won
    points = 0
    # Run through trials
    for trial in trials:
        ###########################################################
        # Check if escape key pressed (or kept pressed since last trial)
        if event.getKeys(['escape']): # exit with Escape
            core.quit()
        # Retrieve picture and answers
        trialnr         = trial['trialnr'] # either retrieve or count yourself
        stimulus        = trial['stimulus']
        condition       = trial['condition']
        valence         = trial['valence']
        reqAction       = trial['reqAction']
        manipulation    = trial['manipulation']
        goValidity      = trial['goValidity']
        nogoValidity    = trial['nogoValidity']
        ISI             = trial['ISI']
        ITI             = trial['ITI']
        print(u'condition is {}'. format(str(condition)))
        print(u'valence is {}'. format(str(valence)))
        print(u'reqAction is {}'. format(str(reqAction)))
        print(u'manipulation is {}'. format(str(manipulation)))
        print(u'goValidity is {}'. format(str(goValidity)))
        print(u'nogoValidity is {}'. format(str(nogoValidity)))
        print(u'ISI is {}'. format(str(ISI)))
        print(u'ITI is {}'. format(str(ITI)))
        ###########################################################
        # Present cue
        if manipulation == 0: # low stakes cue
            stakes = 10 # for points/ outcome later
        elif manipulation == 1: # high stakes cue
            stakes = 50 # for points/ outcome later
            circle_image.draw() # draw circle first (set as image early during initialization)
        else: # 
            error('invalid manipulation value')
        cue_image.setImage('TaskStimuli/{}.jpg'.format(stimulus)) # retrieve cue 
        cue_image.draw() # then draw cue on top
        win.flip()
        ###########################################################
        # Register response
        rtClock = core.Clock()
        keyPress = event.waitKeys(keyList = ['space','escape'], maxWait = 1.3, timeStamped = rtClock)
        time_dif = rtClock.getTime()
        if keyPress:
            if keyPress[0][0] == 'escape': # if abort
                core.quit()
            elif keyPress[0][0] == 'space': # if Go
                response = 1
                RT = time_dif
                print(u'Go!!!!!!!!!')
                print(u'RT is {}'.format(RT))
                core.wait(1.3 - RT) # wait till end of 1.3 sec
            else: 
                print(u'Error in response coding')
                core.quit()
        else: # if NoGo
            response = 0
            RT = 'NA'
        win.flip()
        ###########################################################
        # Check accuracy
        if reqAction == response:
            ACC = 1
        else:
            ACC = 0
        print(u'response is {}'.format(str(response)))
        print(u'ACC is {}'.format(str(ACC)))
        ###########################################################
        # Inter-stimulus interval (ISI) between cue and outcome: 1300 - 1700 ms
        center_fix.setText('+')
        center_fix.draw()
        win.flip()
        core.wait(ISI)
        win.flip()
        ###########################################################
        # Determine validity:
        if response == 1:
            validity = goValidity
        else:
            validity = nogoValidity
        # Select outcome based on cue validity
        if valence == 1 and ACC == 1 and validity == 1:
            outcome = 1;
        elif valence == 1 and ACC == 1 and validity == 0:
            outcome = 0;
        elif valence == 1 and ACC == 0 and validity == 1:
            outcome = 0;
        elif valence == 1 and ACC == 0 and validity == 0:
            outcome = 1;
        elif valence == 0 and ACC == 1 and validity == 1:
            outcome = 0;
        elif valence == 0 and ACC == 1 and validity == 0:
            outcome = -1;
        elif valence == 0 and ACC == 0 and validity == 1:
            outcome = -1;
        elif valence == 0 and ACC == 0 and validity == 0:
            outcome = 0;
        else: 
            outputtext = (u'Error: Cannot determine outcome given valence = {}, ACC = {}, validity = {}').format(str(valence),str(ACC),str(validity))
            print(outputtext)
            core.quit()
        ###########################################################
        # Update points:
        gain = outcome * stakes
        print(u'Gained {} points'.format(str(gain)))
        points += gain
        print(u'Current block score: {} points'.format(str(points)))
        ###########################################################
        # Present outcome 
        if outcome == 1:
            cue_image.setImage('TaskStimuli/reward_{}.jpg'.format(stakes))
        elif outcome == -1:
            cue_image.setImage('TaskStimuli/punishment_{}.jpg'.format(stakes))
        else:
            cue_image.setImage('TaskStimuli/neutral.jpg')
        cue_image.draw()
        win.flip()
        core.wait(0.7) # Present outcome for 700 ms
        win.flip()
        ###########################################################
        # Write results to outcome file (but only of actual trials, not practice trials)
        if block == 'Test':
            file_trials = open('DataOutput/MGNGStakes_{:0>2d}_{:0>2d}.csv'.format(subject,outputAddon),'a') # open output file
            file_trials.write("{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}\n".format(subject,age,gender,hand,part,trialnr,stimulus,condition,reqAction,valence,manipulation,validity,response,ACC,RT,outcome))
            file_trials.close()
        ###########################################################
        # Intertrial interval (ITI): 1800-2300 ms
        center_fix.setText('+')
        center_fix.draw()
        win.flip()
        core.wait(ITI) # wait for ITI (in seconds) so pupil returns to baseline levels
        win.flip()
        ###########################################################
        # End of trials
    # Calculate global performance statistics
    return(points)

###########################################################################
# Run block including instructions and practice stimuli
###########################################################################

def function_block(block = None, part = None):
    # Load stimuli with trialhandler
    if block == 'Pract':
        stimulusList = data.importConditions("TrialHandler/stimuluslist_pract_Part{}.csv".format(part)) # Load test stimuli   
    elif block == 'Test':
        stimulusList = data.importConditions("TrialHandler/stimuluslist_test_stimNum_{}_sub_{:0>3d}_Part{}.csv".format(nTrials,subject, part)) # Load test stimuli   
    else:
        outputtext = (u'Error: Invalid settings for coming block.')
        print(outputtext)
        core.quit()        
    # Prepare trials:    
    trials = data.TrialHandler(stimulusList, nReps = 1, method = "sequential")
    # Do task:
    (points) = function_run_trials(trials = trials, block = block, part = part)
    return(points)
###################################################################################
# Main function (allows for running the program if executed from the file itself)
###################################################################################

def function_main():
    # Display instructions at beginning
    if isPractice == True:
        # General instructions:
        function_display_instructions('Instructions/MGNGStakes_Instructions_General_01.jpg')
        function_display_instructions('Instructions/MGNGStakes_Instructions_General_02.jpg')
        function_display_instructions('Instructions/MGNGStakes_Instructions_General_03.jpg')
        function_display_instructions('Instructions/MGNGStakes_Instructions_General_04.jpg')
        function_display_instructions('Instructions/MGNGStakes_Instructions_General_05.jpg')
        function_display_instructions('Instructions/MGNGStakes_Instructions_General_06.jpg')
        function_display_instructions('Instructions/MGNGStakes_Instructions_General_07.jpg')
        # Practice rounds:
        function_display_instructions('Instructions/MGNGStakes_Instructions_Pract_G2W01.jpg')
        function_block(block = 'Pract', part = 1) # 4 Go2Win trials
        function_display_instructions('Instructions/MGNGStakes_Instructions_Pract_G2W02.jpg')
        function_display_instructions('Instructions/MGNGStakes_Instructions_Pract_NG2W01.jpg')
        function_block(block = 'Pract', part = 2) # 4 NoGo2Win trials
        function_display_instructions('Instructions/MGNGStakes_Instructions_Pract_NG2W02.jpg')
        function_display_instructions('Instructions/MGNGStakes_Instructions_Pract_G2A01.jpg')
        function_block(block = 'Pract', part = 3) # 4 Go2Avoid trials
        function_display_instructions('Instructions/MGNGStakes_Instructions_Pract_G2A02.jpg')
        function_display_instructions('Instructions/MGNGStakes_Instructions_Pract_NG2A01.jpg')
        function_block(block = 'Pract', part = 4) # 4 NoGo2Avoid trials
        function_display_instructions('Instructions/MGNGStakes_Instructions_Pract_NG2A02.jpg')
#   # Test Rounds:
    allPoints = []; # empty list that will contain points of each block
    function_display_instructions('Instructions/MGNGStakes_Instructions_BeforeTest_01.jpg')
    function_display_instructions('Instructions/MGNGStakes_Instructions_BeforeTest_02.jpg')
    if isBlock1 == True:
        # Block 1:
        function_display_instructions('Instructions/MGNGStakes_Instructions_BeforeTest_Block01.jpg')
        function_display_instructions('Instructions/MGNGStakes_Instructions_Remember.jpg')
        (points) = function_block(block = 'Test', part = 1)
        allPoints.append(points)
        print(u'>>> Current total score after block 1: {} points'.format(str(points)))
        function_display_instructions('Instructions/MGNGStakes_Instructions_Break_Block01.jpg')
    if isBlock2 == True:
        # Block 2:
        function_display_instructions('Instructions/MGNGStakes_Instructions_BeforeTest_Block02.jpg')
        function_display_instructions('Instructions/MGNGStakes_Instructions_Remember.jpg')
        (points) = function_block(block = 'Test', part = 2)
        allPoints.append(points)
        print(u'>>> Current total score after block 2: {} points'.format(str(points)))
        function_display_instructions('Instructions/MGNGStakes_Instructions_Break_Block02.jpg')
    if isBlock3 == True:
        # Block 3:
        function_display_instructions('Instructions/MGNGStakes_Instructions_BeforeTest_Block03.jpg')
        function_display_instructions('Instructions/MGNGStakes_Instructions_Remember.jpg')
        (points) = function_block(block = 'Test', part = 3)
        allPoints.append(points)
        print(u'>>> Current total score after block 3: {} points'.format(str(points)))
        function_display_instructions('Instructions/MGNGStakes_Instructions_Break_Block03.jpg')
        # Block 4:
    if isBlock4 == True:
        function_display_instructions('Instructions/MGNGStakes_Instructions_BeforeTest_Block04.jpg')
        function_display_instructions('Instructions/MGNGStakes_Instructions_Remember.jpg')
        (points) = function_block(block = 'Test', part = 4)
        allPoints.append(points)
        print(u'>>> Current total score after block 4: {} points'.format(str(points)))
    # Compute sum of all points in all blocks:
    totalPoints = sum(allPoints)
    # Displays scores
    # First print (for security in case participants advance task too fast)
    print(u'The task is over now. You have gained {} points.\nPlease contact the experimenter now.').format(str(totalPoints))
    # Then on screen
    center_text.setPos((0,.2)) # first line of text
    center_text.setText(u'The task is over now. You have gained {} points.'.format(str(totalPoints)))
    center_text.draw()
    center_text.setPos((0,0)) # second line of text
    if totalPoints >= nTrials*3: # 960 points with 20 reps; equivalent to 66.67% accuracy
        center_text.setText(u'This gives you an extra reward!')
    else:
        center_text.setText(u'Unfortunately, that is not enough for an extra reward.')
    center_text.draw()
    center_text.setPos((0,-.2)) # third line of text
    center_text.setText(u'Please contact the experimenter now.')
    center_text.draw()
    win.flip()
    waitforbutton(waitTime = 2) # prevent participants from advancing task prematurely 
    win.flip() # clean screen 
    function_display_instructions('Instructions/MGNGStakes_Instructions_End01.jpg')
###################################################################################################################################################################################
###################################################################################################################################################################################
###################################################################################################################################################################################

###########################################################################
# Get initial participant info
(subject, age, gender, hand, isPractice, isBlock1, isBlock2, isBlock3, isBlock4, nTrials) = function_gui_info()
# Create output file addon
outputAddon = random.randrange(0,100)

###########################################################################
# Create global variables

# Window:
# Mind: color in PsychoPy is between -1 and 1, so use rgb/355*2-1 for conversion
win=visual.Window(fullscr=True, color=[.30,.30,.30], winType='pygame', units = 'norm') # check gray background color with PowerPoint; is 166/255=0.6509804 in original stimuli
win.setMouseVisible(False)
# Central text (for fixation cross and final points feedback):
center_fix  = visual.TextStim(win=win,text="", pos = (0,0), color=[0.03,0.03,0.03], wrapWidth = 0.8, alignHoriz = 'center', height = .2) # Test out color: .80? .30
center_text = visual.TextStim(win=win,text="", pos = (0,0), color=[1,1,1], wrapWidth = 0.8, alignHoriz = 'center', height = .1) # Test out color: white enough? used for final feedback
# Images:
instr_image = visual.ImageStim(win=win,image='Instructions/MGNGStakes_Instructions_General_01.jpg', units = 'norm')
cue_image = visual.ImageStim(win=win,image='TaskStimuli/A1.jpg', pos = (0,0), size = (300, 300), units = 'pix') # original 720, 720 pixels
circle_image = visual.ImageStim(win=win,image='TaskStimuli/circle_dondersred.png', pos = (0,0), size = (500, 500), units = 'pix') # original 720, 720 pixels
 
# Define trial output file
file_trials = open('DataOutput/MGNGStakes_{:0>2d}_{:0>2d}.csv'.format(subject,outputAddon),'w') # write headers of file: 
file_trials.write("{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}\n".format('Subject','Age','Gender','Hand','Block','Trialnr','Stimulus','Condition','ReqAction','Valence','Manipulation','Validity','Response','ACC','RT','Outcome'))
file_trials.close()

###########################################################################
# Run experiment
###########################################################################

function_main()

###########################################################################
# End experiment
###########################################################################
core.quit()
# END of program