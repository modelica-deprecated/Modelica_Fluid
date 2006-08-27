package ControlledTankSystem 
  "Tank system with controller, start/stop/shut operation and diagram animation" 
  model ControlledTanks 
    "Demonstrating the controller of a tank filling/emptying system" 
    extends Modelica.Icons.Example;
    package Medium = Modelica.Media.Water.ConstantPropertyLiquidWater;
    
    Modelica_Fluid.Examples.ControlledTankSystem.Utilities.TankController 
      tankController(
      waitTime=50,
      minLevel=0.001,
      maxLevel=0.9*tank1.levelMax) 
      annotation (extent=[-60,-20; -20,20]);
    Modelica_Fluid.Examples.ControlledTankSystem.Utilities.RadioButton start(
                                                           reset={stop.on,shut.on},
        buttonTimeTable={20,280}) 
      annotation (extent=[-100,20; -80,40]);
    Modelica_Fluid.Examples.ControlledTankSystem.Utilities.RadioButton stop(
                                                          reset={start.on,shut.on},
        buttonTimeTable={220,650}) 
      annotation (extent=[-100,-10; -80,10]);
    Modelica_Fluid.Examples.ControlledTankSystem.Utilities.RadioButton shut(
                                                          reset={start.on,stop.on},
        buttonTimeTable={700}) 
      annotation (extent=[-100,-40; -80,-20]);
    annotation (
      Diagram,
      Coordsys(grid=[1,1],      component=[20, 20],
        scale=0),
      experiment(StopTime=900),
      experimentSetupOutput,
      Documentation(info="<html>
<p>
With this example, the controller of a tank filling/emptying system
is demonstrated. 
</p>

<p>
The basic operation is to fill and empty the two tanks:
</p>
<ol>
<li> Valve 1 is opened and tank 1 is filled.</li>
<li> When tank 1 reaches its fill level limit, 
     valve 1 is closed. </li>
<li> After a waiting time, valve 2 is
     opened and the fluid flows from tank 1 into tank 2.</li>
<li> When tank 1 reaches its minimum level, valve 2 is closed. </li>
<li> After a waiting time, valve 3 is opened and
     the fluid flows out of tank 2</li>
<li> When tank 2 reaches its minimum level, valve 3 is closed</liI>
</ol>
<p>
The above \"normal\" process can be influenced by three 
buttons:
</p>
<ul>
<li> Button <b>start</b> starts the above process.
     When this button is pressed after a \"stop\" or
     \"shut\" operation, the process operation continues.
     </li>.
<li> Button <b>stop</b> stops the above process by
     closing all valves. Then, the controller waits for
     further input (either \"start\" or \"shut\" operation).</li>
<li> Button <b>shut</b> is used to shutdown the process, 
     by emptying at once both tanks by opening valve 2 and
     valve 3. When this is achieved,
     the process goes back to its start configuration
     where all 3 valves are closed.
     Clicking on \"start\", restarts the process.</li>
</ul> 
     
<p>
The demo-run uses the following button presses:
</p>
 
<ul>
<li> Button <b>start</b> pressed at 20 s.</li>
<li> Button <b>stop</b> pressed at 220 s </li>
<li> Button <b>start</b> pressed at 280 s </li>
<li> Button <b>stop</b> pressed at 650 s </li>
<li> Button <b>shut</b> pressed at 700 s </li>
<li> Simulate for 900 s</li>
</ul>

<p>
This example is based on
</p>

<dl>
<dt>Dressler I. (2004):</dt>
<dd> <b>Code Generation From JGrafchart to Modelica</b>.
     Master thesis, supervisor: Karl-Erik Arzen,
     Department of Automatic Control, Lund Institute of Technology,
     Lund, Sweden, March 30, 2004<br>&nbsp;</dd>
</dl>
 
</html>"),
      Commands(file=
            "../Scripts/Examples/ControlledTanks/plot level and ports.m_flow.mos" 
          "plot level and ports.m_flow"));
    ControlValves.ValveDiscrete valve1(                             redeclare 
        package Medium = Medium, Kv=4e-4) 
      annotation (extent=[10,40; 30,60], rotation=90);
    Volumes.Tank tank1(
      level_start=0.05,
      redeclare package Medium = Medium,
      nTopPorts=1,
      portsData={Modelica_Fluid.Volumes.BaseClasses.TankPortData(diameter=0.2,
          portLevel=0)},
      area=6,
      levelMax=4)             annotation (extent=[0,-10; 40,30]);
    Modelica.Blocks.Sources.RealExpression level1(y=tank1.level) 
      annotation (extent=[-90,-60; -55,-40]);
    ControlValves.ValveDiscrete valve2(                redeclare package Medium
        = Medium, Kv=100) 
      annotation (extent=[10,-40; 30,-20], rotation=90);
    ControlValves.ValveDiscrete valve3(                redeclare package Medium
        = Medium, Kv=10) 
      annotation (extent=[70,-80; 90,-60], rotation=90);
    Volumes.Tank tank2(
      level_start=0.05,
      redeclare package Medium = Medium,
      levelMax=5,
      area=6,
      nTopPorts=1,
      portsData={Modelica_Fluid.Volumes.BaseClasses.TankPortData(diameter=0.2,
          portLevel=0)})      annotation (extent=[60,-50; 100,-10]);
    Modelica_Fluid.Sources.FixedBoundary_pTX ambient1(
                                                     redeclare package Medium 
        = Medium, 
      p=ambient.default_p_ambient, 
      T=ambient.default_T_ambient) 
      annotation (extent=[10,-100; 30,-80], rotation=0);
    Modelica.Blocks.Sources.RealExpression level2(y=tank2.level) 
      annotation (extent=[-70,-80; -33,-60]);
    Modelica_Fluid.Sources.FixedBoundary_pTX source(
                                                   redeclare package Medium = 
          Medium, p=2.5e6, 
      T=ambient.default_T_ambient) 
      annotation (extent=[10,70; 30,90],    rotation=-90);
    inner Ambient ambient annotation (extent=[-90,70; -70,90]);
  equation 
    connect(shut.on, tankController.shut) annotation (points=[-79,-30; -72,-30;
          -72,-12; -62,-12],    style(color=5, rgbcolor={255,0,255}));
    connect(stop.on, tankController.stop) annotation (points=[-79,0; -62,0],
        style(color=5, rgbcolor={255,0,255}));
    connect(start.on, tankController.start) annotation (points=[-79,30; -70,30;
          -70,12; -62,12],       style(color=5, rgbcolor={255,0,255}));
    connect(tankController.valve1, valve1.open) annotation (points=[-19,12; -10,
          12; -10,50; 12,50], style(color=5, rgbcolor={255,0,255}));
    connect(level1.y, tankController.level1) annotation (points=[-53.25,-50;
          -52,-50; -52,-22],
                    style(color=74, rgbcolor={0,0,127}));
    connect(tankController.valve2, valve2.open) annotation (points=[-19,0; -5,0;
          -5,-30; 12,-30], style(color=5, rgbcolor={255,0,255}));
    connect(tankController.valve3, valve3.open) annotation (points=[-19,-12;
          -10,-12; -10,-70; 72,-70],
                                 style(color=5, rgbcolor={255,0,255}));
    connect(level2.y, tankController.level2) annotation (points=[-31.15,-70;
          -28,-70; -28,-22],
                    style(color=74, rgbcolor={0,0,127}));
    connect(source.port, valve1.port_b) 
      annotation (points=[20,70; 20,60], style(color=69, rgbcolor={0,127,255}));
    connect(valve1.port_a, tank1.topPorts[1]) 
      annotation (points=[20,40; 20,30], style(color=69, rgbcolor={0,127,255}));
    connect(tank1.ports[1], valve2.port_b) annotation (points=[20,-10; 20,-20],
        style(color=69, rgbcolor={0,127,255}));
    connect(valve2.port_a, tank2.topPorts[1]) annotation (points=[20,-40; 20,
          -50; 50,-50; 50,2; 80,2; 80,-10],style(color=69, rgbcolor={0,127,255}));
    connect(tank2.ports[1], valve3.port_b) annotation (points=[80,-50; 80,-60], style(
          color=69, rgbcolor={0,127,255}));
    connect(valve3.port_a, ambient1.port) annotation (points=[80,-80; 80,-90;
          30,-90], style(color=69, rgbcolor={0,127,255}));
    
  end ControlledTanks;
  
  package Utilities 
    model TankController "Controller for tank system" 
      import SI = Modelica.SIunits;
      extends Modelica.StateGraph.Interfaces.PartialStateGraphIcon;
      parameter SI.Height maxLevel "Fill level of tank 1";
      parameter SI.Height minLevel "Lowest level of tank 1 and 2";
      parameter SI.Time waitTime "Wait time, between operations";
      
      Modelica.StateGraph.InitialStep s1(nIn=2) 
                     annotation (extent=[-72,30; -52,50]);
      Modelica_Fluid.Examples.ControlledTankSystem.Utilities.NormalOperation 
        normal(
        maxLevel=maxLevel,
        minLevel=minLevel,
        waitTime=waitTime) 
        annotation (extent=[-20,20; 20,60]);
      Modelica.StateGraph.Transition T1(condition=start) 
                                     annotation (extent=[-50,50; -30,30]);
      Modelica.StateGraph.Transition T2(condition=level2 < minLevel) 
        annotation (extent=[27,50; 47,30]);
      Modelica.StateGraph.Transition T3(condition=stop) 
        annotation (extent=[-33,-11; -13,9],     rotation=-90);
      Modelica.StateGraph.Step s2(nOut=2) 
              annotation (extent=[-50,-60; -30,-40]);
      Modelica.StateGraph.Transition T4(condition=start) 
        annotation (extent=[0,-10; 20,10],    rotation=90);
      Modelica.StateGraph.Transition T5(condition=shut) 
                                    annotation (extent=[-6,-60; 14,-40]);
      Modelica.StateGraph.Step emptyTanks 
                      annotation (extent=[20,-60; 40,-40]);
      Modelica.StateGraph.Transition T6(condition=level1 < minLevel and level2
             < minLevel) 
        annotation (extent=[45,-60; 65,-40]);
      Modelica.Blocks.Interfaces.BooleanInput start 
        annotation (extent=[-120, 50; -100, 70]);
      Modelica.Blocks.Interfaces.BooleanInput stop 
        annotation (extent=[-120, -10; -100, 10]);
      Modelica.Blocks.Interfaces.BooleanInput shut 
        annotation (extent=[-120, -70; -100, -50]);
      Modelica.Blocks.Interfaces.RealInput level1 
        annotation (extent=[-70, -120; -50, -100], rotation=90);
      Modelica.Blocks.Interfaces.RealInput level2 
        annotation (extent=[50, -120; 70, -100], rotation=90);
      Modelica.Blocks.Interfaces.BooleanOutput valve1 
        annotation (extent=[100, 55; 110, 65]);
      Modelica.Blocks.Interfaces.BooleanOutput valve2 
        annotation (extent=[100, -5; 110, 5]);
      Modelica.Blocks.Interfaces.BooleanOutput valve3 
        annotation (extent=[100, -65; 110, -55]);
      Modelica.Blocks.Sources.BooleanExpression setValve1(y=normal.fillTank1.
            active) 
        annotation (extent=[20,73; 80,92]);
      Modelica.Blocks.Sources.BooleanExpression setValve2(y=normal.fillTank2.
            active or emptyTanks.active and level1 > minLevel) 
        annotation (extent=[-40,-85; 80,-64]);
      Modelica.Blocks.Sources.BooleanExpression setValve3(y=normal.emptyTank2.
            active or emptyTanks.active and level2 > minLevel) 
        annotation (extent=[-40,-103; 80,-83]);
    equation 
      
      annotation (structurallyIncomplete,
        Diagram(Rectangle(extent=[-100,100; 100,-100],   style(color=0,
                rgbcolor={0,0,0}))),
        Icon(
          Text(
            extent=[-100,68; -12,54],
            style(
              color=0,
              fillColor=0,
              fillPattern=1),
            string="start"),
          Text(
            extent=[-99,6; -14,-9],
            style(
              color=0,
              fillColor=0,
              fillPattern=1),
            string="stop"),
          Text(
            extent=[-99,-54; -14,-69],
            style(
              color=0,
              fillColor=0,
              fillPattern=1),
            string="shut"),
          Text(
            extent=[-94,-82; -9,-96],
            style(
              color=0,
              fillColor=0,
              fillPattern=1),
            string="level1"),
          Text(
            extent=[11,-83; 96,-98],
            style(
              color=0,
              fillColor=0,
              fillPattern=1),
            string="level2"),
          Text(
            extent=[10,68; 99,54],
            style(
              color=0,
              fillColor=0,
              fillPattern=1),
            string="valve1"),
          Text(
            extent=[7,10; 101,-5],
            style(
              color=0,
              fillColor=0,
              fillPattern=1),
            string="valve2"),
          Text(
            extent=[2,-51; 102,-67],
            style(
              color=0,
              fillColor=0,
              fillPattern=1),
            string="valve3")),
        Coordsys(grid=[1,1],  component=[20, 20],
          scale=0));
      connect(s1.outPort[1], T1.inPort) 
                                     annotation (points=[-51.5,40; -44,40],
          style(
          color=0,
          fillColor=0,
          fillPattern=1));
      connect(T1.outPort, normal.inPort)      annotation (points=[-38.5,40;
            -21.3333,40],
                  style(
          color=0,
          fillColor=0,
          fillPattern=1));
      connect(normal.outPort, T2.inPort)      annotation (points=[20.6667,40;
            33,40],
                  style(
          color=0,
          fillColor=0,
          fillPattern=1));
      connect(T5.outPort, emptyTanks.inPort[1]) 
                                             annotation (points=[5.5,-50; 19,
            -50],  style(
          color=0,
          fillColor=0,
          fillPattern=1));
      connect(emptyTanks.outPort[1], T6.inPort) 
                                             annotation (points=[40.5,-50; 51,
            -50],  style(
          color=0,
          fillColor=0,
          fillPattern=1));
      connect(setValve1.y, valve1) 
        annotation (points=[83,82.5; 90,82.5; 90,60; 105,60], style(color=5));
      connect(setValve2.y, valve2) 
        annotation (points=[86,-74.5; 90,-74.5; 90,0; 105,0], style(color=5));
      connect(setValve3.y, valve3) annotation (points=[86,-93; 95,-93; 95,-60;
            105,-60],   style(color=5));
      connect(normal.suspend[1], T3.inPort)   annotation (points=[-10,19.3333;
            -10,12; -23,12; -23,3],  style(color=0, rgbcolor={0,0,0}));
      connect(T3.outPort, s2.inPort[1]) 
                                     annotation (points=[-23,-2.5; -23,-20; -60,
            -20; -60,-50; -51,-50],        style(color=0, rgbcolor={0,0,0}));
      connect(level1, normal.level1)      annotation(points=[-60,-110; -60,-80;
            -80,-80; -80,20; -30,20; -30,24; -22.6667,24],
                                                      style(color=3, rgbcolor={
              0,0,255}));
      connect(s2.outPort[1], T5.inPort) annotation(points=[-29.5,-49.75; -30,
            -49.75; -30,-50; 0,-50], style(color=0, rgbcolor={0,0,0}));
      connect(s2.outPort[2], T4.inPort) annotation(points=[-29.5,-50.25; -29,-50;
            -8,-50; -8,-25; 10,-25; 10,-4], style(color=0, rgbcolor={0,0,0}));
      connect(T2.outPort, s1.inPort[1]) annotation(points=[38.5,40; 70,40; 70,
            70; -80,70; -80,40; -73,40; -73,40.5],
                                             style(color=0, rgbcolor={0,0,0}));
      connect(T6.outPort, s1.inPort[2]) annotation(points=[56.5,-50; 70,-50; 70,
            70; -80,70; -80,40; -74,40; -73,39.5],
                                             style(color=0, rgbcolor={0,0,0}));
      connect(T4.outPort, normal.resume[1])      annotation (points=[10,1.5; 10,
            10; 10.5,10; 10.5,18.6667; 10,18.6667],  style(color=0, rgbcolor={0,0,0}));
    end TankController;
    
    model NormalOperation 
      "Normal operation of tank system (button start pressed)" 
      import SI = Modelica.SIunits;
      extends Modelica.StateGraph.PartialCompositeStep;
      parameter SI.Height maxLevel "Fill level of tank 1";
      parameter SI.Height minLevel "Lowest level of tank 1 and 2";
      parameter SI.Time waitTime "Wait time between operations";
      
      Modelica.Blocks.Interfaces.RealInput level1 
        annotation (extent=[-190,-140; -150,-100]);
      annotation (Diagram, Documentation(info="<html>
 
</html>"));
      Modelica.StateGraph.Step fillTank1 
                     annotation (extent=[-140,-10; -120,10]);
      Modelica.StateGraph.Transition T1(condition=level1 > maxLevel) 
        annotation (extent=[-110,-10; -90,10]);
      Modelica.StateGraph.Step fillTank2 
                     annotation (extent=[-10,-10; 10,10]);
      Modelica.StateGraph.Transition T3(condition=level1 < minLevel) 
        annotation (extent=[20,-10; 40,10]);
      Modelica.StateGraph.Step emptyTank2 
                      annotation (extent=[120,-10; 140,10]);
      Modelica.StateGraph.Step wait1 
                 annotation (extent=[-80,-10; -60,10]);
      Modelica.StateGraph.Transition T2(enableTimer=true, waitTime=waitTime) 
        annotation (extent=[-50,-10; -30,10]);
      Modelica.StateGraph.Step wait2 
                 annotation (extent=[54,-10; 74,10]);
      Modelica.StateGraph.Transition T4(enableTimer=true, waitTime=waitTime) 
        annotation (extent=[82,-10; 102,10]);
    equation 
      connect(fillTank1.inPort[1], inPort) 
                                        annotation (points=[-141,0; -160,0],
          style(
          color=0,
          fillColor=0,
          fillPattern=1));
      connect(fillTank1.outPort[1], T1.inPort) 
                                            annotation (points=[-119.5,0; -104,
            0],  style(
          color=0,
          fillColor=0,
          fillPattern=1));
      connect(fillTank2.outPort[1], T3.inPort) 
                                            annotation (points=[10.5,0; 26,0],
           style(
          color=0,
          fillColor=0,
          fillPattern=1));
      connect(emptyTank2.outPort[1], outPort) 
                                           annotation (points=[140.5,0; 155,0],
           style(
          color=0,
          fillColor=0,
          fillPattern=1));
      connect(wait1.outPort[1], T2.inPort) 
                                        annotation (points=[-59.5,0; -44,0],
          style(color=0, rgbcolor={0,0,0}));
      connect(T2.outPort, fillTank2.inPort[1]) 
                                            annotation (points=[-38.5,0; -11,0],
                style(color=0, rgbcolor={0,0,0}));
      connect(T1.outPort, wait1.inPort[1]) 
                                        annotation (points=[-98.5,0; -81,0],
           style(color=0, rgbcolor={0,0,0}));
      connect(wait2.outPort[1], T4.inPort) 
                                        annotation (points=[74.5,0; 88,0],
          style(color=0, rgbcolor={0,0,0}));
      connect(T3.outPort, wait2.inPort[1]) 
        annotation (points=[31.5,0; 53,0],   style(color=0, rgbcolor={0,0,0}));
      connect(T4.outPort,emptyTank2.inPort[1]) 
                                             annotation (points=[93.5,0; 119,0],
                 style(color=0, rgbcolor={0,0,0}));
    end NormalOperation;
    
    block RadioButton 
      "Button that sets its output to true when pressed and is reset when an element of 'reset' becomes true" 
      
      parameter Modelica.SIunits.Time buttonTimeTable[:] 
        "Time instants where button is pressend";
      input Boolean reset[:]={false} 
        "Reset button to false, if an element of reset becomes true" 
        annotation (Dialog(group="Time varying expressions"));
      
      annotation (Icon(
          onClick=setVariable(on, true),
          Rectangle(extent=[-100, -100; 100, 100], style(
              color=10,
              thickness=2,
              fillColor=DynamicSelect(8, if on > 0.5 then 2 else 8),
              fillPattern=DynamicSelect(11, if on > 0.5 then 12 else 11))),
          Text(
            extent=[-80, -40; 80, 40],
            style(color=0),
            string="%name")), Diagram,
          Documentation(info="<html>
  
</html>"));
      Modelica.Blocks.Interfaces.BooleanOutput on 
        annotation (extent=[100, -10; 120, 10], style(color=0));
    protected 
      Modelica.Blocks.Sources.BooleanTable table(table=buttonTimeTable);
    algorithm 
      when pre(reset) then
         on := false;
      end when;
      
      when change(table.y) then
         on := true;
      end when;
    end RadioButton;
  end Utilities;
end ControlledTankSystem;
