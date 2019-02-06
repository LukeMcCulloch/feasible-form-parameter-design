#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 27 22:37:42 2017

@author: luke

http://chatterbot.readthedocs.io/en/latest/training.html


You can also specify file paths to corpus files or directories of corpus files when calling the train method.

chatterbot.train(
    "./data/greetings_corpus/custom.corpus.json",
    "./data/my_corpus/"
)
"""

# -*- coding: utf-8 -*-
from chatterbot import ChatBot
from chatterbot.trainers import ChatterBotCorpusTrainer


chatterbot = ChatBot(
    "Math & Time Bot",
    logic_adapters=[
        "chatterbot.logic.MathematicalEvaluation",
        "chatterbot.logic.TimeLogicAdapter"
    ],
    input_adapter="chatterbot.input.VariableInputTypeAdapter",
    output_adapter="chatterbot.output.OutputAdapter"
)




# Print an example of getting one math based response
response = chatterbot.get_response("What is 4 + 9?")
print(response)

# Print an example of getting one time based response
response = chatterbot.get_response("What time is it?")
print(response)


# Print an example of getting one time based response
response = chatterbot.get_response("I am going to train you...")
print(response)

chatterbot.set_trainer(ChatterBotCorpusTrainer)

#chatterbot.train(
#    "chatterbot.corpus.english"
#)


chatterbot.train(
    "./bot_project_knowledge/basic_overview/ProjectInfo.yml"
)



# Print an example of getting one time based response
response = chatterbot.get_response("What is Luke's PhD project?")
print(response)

"""
from chatterbot import ChatBot
from chatterbot.trainers import ListTrainer


chatterbot = ChatBot("Ron Obvious")


conversation = [
    "Hello",
    "Hi there!",
    "How are you doing?",
    "I'm doing great.",
    "That is good to hear",
    "Thank you.",
    "You're welcome."
]

chatterbot.set_trainer(ListTrainer)
chatterbot.train(conversation)

response = chatterbot.get_response("Good morning!")
print(response)
"""