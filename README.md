# Flight-path-reconstruction-state-and-fault-prediction-model
This project focused on the computer’s ability to reconstruct, predict and understand the flight path and state of a drone using biased, incomplete and unclean sensor data.  End result was to achieve the automatic prediction of the cause of a drone accident.

Given that loss of control is a vital part of this study, it has to be fully defined. A general definition of loss of control is the significant departure of the aircraft from the flight envelope. Faults investigated were:

1. Sensor blackout - Sensor outages occur when a sensor stops recording data, this typically happens to GPS
sensors.
2. Jammed control surface - occurs the control surface is jammed at any angle of deflection less than the maximum. 
3. Control surface hard-over - occurs when the control surface is jammed at the maximum angle of deflection.
4. Control surface float - occurs either hydraulics or severing or loosening of cables may result in free oscillations of control surfaces, also known as control surface float

#An overview of one of the Kalman Filters is shown below:

![image](https://user-images.githubusercontent.com/55096537/111266265-f31a3580-8632-11eb-9a94-47d0cfbb7d57.png)

