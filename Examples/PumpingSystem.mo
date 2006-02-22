model PumpingSystem "Model of a pumping system for drinking water" 
  extends Modelica.Icons.Example;
  model HeatExchanger "Double pipe heat exchanger with outer wall neglected" 
    
    //General
    parameter Integer n(min=1) "Spatial segmentation";
    replaceable package Medium_1 = Modelica.Media.Water.StandardWater extends 
      Modelica.Media.Interfaces.PartialMedium "Inner pipe medium" 
                                                      annotation(choicesAllMatching);
    replaceable package Medium_2 = Modelica.Media.Water.StandardWater extends 
      Modelica.Media.Interfaces.PartialMedium "Outer pipe medium" 
                                                      annotation(choicesAllMatching);
    parameter SI.Length di_1(min=0) "Inner diameter of inner pipe"     annotation(Dialog(tab="General", group="Dimensions"));
    parameter SI.Length da_1(min=0) "Outer diameter of inner pipe"     annotation(Dialog(tab="General", group="Dimensions"));
    parameter SI.Length di_2(min=0) "Inner diameter of outer pipe"     annotation(Dialog(tab="General", group="Dimensions"));
    parameter SI.Length length(min=0) "Length of both pipes" annotation(Dialog(tab="General", group="Dimensions"));
    
    //Wall
    parameter SI.Density d_wall "Density of wall material" annotation(Dialog(tab="General", group="Constant material properties"));
    parameter SI.SpecificHeatCapacity c_wall 
      "Specific heat capacity of wall material"                                        annotation(Dialog(tab="General", group="Constant material properties"));
    final parameter SI.Mass m_wall=sum(pipe_1.wall.m) "Wall mass";
    parameter Boolean initWall_steadyState=false 
      "= true, Wall initialization in steady state"                                      annotation(Dialog(tab="Initialization", group="Wall"));
    parameter SI.Temperature T_start_wall "Start value of wall temperature" annotation(Dialog(tab="Initialization", group="Wall"));
    
    //Initialization pipe 1
    parameter Modelica_Fluid.Types.Init.Temp initOption_1 
      "Initialization option" 
      annotation(Evaluate=true, Dialog(tab = "Initialization", group = "Inner pipe"));
    parameter Boolean use_T_start_1=true 
      "Use T_start if true, otherwise h_start" 
      annotation(Evaluate=true, Dialog(tab = "Initialization", group = "Inner pipe"));
    parameter Medium_1.AbsolutePressure p_a_start1=Medium_1.p_default 
      "Start value of pressure" 
      annotation(Dialog(tab = "Initialization", group = "Inner pipe"));
    parameter Medium_1.AbsolutePressure p_b_start1=Medium_1.p_default 
      "Start value of pressure" 
      annotation(Dialog(tab = "Initialization", group = "Inner pipe"));
    parameter Medium_1.Temperature T_start_1=if use_T_start_1 then Medium_1.T_default else 
        Medium_1.temperature_phX((p_a_start1+p_b_start1)/2, h_start_1, X_start_1) 
      "Start value of temperature" 
      annotation(Evaluate=true, Dialog(tab = "Initialization", group = "Inner pipe", enable = use_T_start_1));
    parameter Medium_1.SpecificEnthalpy h_start_1=if use_T_start_1 then 
        Medium_1.specificEnthalpy_pTX((p_a_start1+p_b_start1)/2, T_start_1, X_start_1[1:Medium_1.nXi]) else Medium_1.h_default 
      "Start value of specific enthalpy" 
      annotation(Evaluate=true, Dialog(tab = "Initialization", group = "Inner pipe", enable = not use_T_start_1));
    parameter Medium_1.MassFraction X_start_1[Medium_1.nX]=Medium_1.X_default 
      "Start value of mass fractions m_i/m" 
      annotation (Dialog(tab="Initialization", group = "Inner pipe", enable=(Medium_1.nXi > 0)));
    parameter Medium_1.MassFlowRate mflow_start_1 
      "Start value of mass flow rate" annotation(Evaluate=true, Dialog(tab = "Initialization", group = "Inner pipe"));
    //Initialization pipe 2
    parameter Modelica_Fluid.Types.Init.Temp initOption_2 
      "Initialization option" 
      annotation(Evaluate=true, Dialog(tab = "Initialization", group = "Outer pipe"));
    parameter Boolean use_T_start_2=true 
      "Use T_start if true, otherwise h_start" 
      annotation(Evaluate=true, Dialog(tab = "Initialization", group = "Outer pipe"));
    parameter Medium_2.AbsolutePressure p_a_start2=Medium_2.p_default 
      "Start value of pressure" 
      annotation(Dialog(tab = "Initialization", group = "Outer pipe"));
    parameter Medium_2.AbsolutePressure p_b_start2=Medium_2.p_default 
      "Start value of pressure" 
      annotation(Dialog(tab = "Initialization", group = "Outer pipe"));
    parameter Medium_2.Temperature T_start_2=if use_T_start_2 then Medium_2.T_default else 
        Medium_2.temperature_phX((p_a_start2+p_b_start2)/2, h_start_2, X_start_2) 
      "Start value of temperature" 
      annotation(Evaluate=true, Dialog(tab = "Initialization", group = "Outer pipe", enable = use_T_start_2));
    parameter Medium_2.SpecificEnthalpy h_start_2=if use_T_start_2 then 
        Medium_2.specificEnthalpy_pTX((p_a_start2+p_b_start2)/2, T_start_2, X_start_2[1:Medium_2.nXi]) else Medium_2.h_default 
      "Start value of specific enthalpy" 
      annotation(Evaluate=true, Dialog(tab = "Initialization", group = "Outer pipe", enable = not use_T_start_2));
    parameter Medium_2.MassFraction X_start_2[Medium_2.nX]=Medium_2.X_default 
      "Start value of mass fractions m_i/m" 
      annotation (Dialog(tab="Initialization", group = "Outer pipe", enable=Medium_2.nXi>0));
    parameter Medium_2.MassFlowRate mflow_start_2 
      "Start value of mass flow rate"    annotation(Evaluate=true, Dialog(tab = "Initialization", group = "Outer pipe"));
    //Advanced
    parameter Boolean lumped_dp = true 
      " = true, lumped pressure drop, reduces number of pressure states to one"
                                annotation(Evaluate=true, Dialog(tab="Advanced", group="Conservation of mass, energy, momentum"));
    parameter Boolean static "= true, use quasistatic mass and energy balances"
                             annotation(Evaluate=true, Dialog(tab="Advanced", group="Conservation of mass, energy, momentum"));
    parameter Boolean kineticTerm 
      " = true, include kinetic term in momentum balance" 
                                  annotation(Evaluate=true, Dialog(tab="Advanced", group="Conservation of mass, energy, momentum"));
    parameter Real K1=1 
      "Enhancement factor for heat transfer area pipe 1(=>parallel tubes)"  annotation(Dialog(tab="General", group="Heat transfer"));
    parameter Real K2=1 
      "Enhancement factor for heat transfer area pipe 2(=>parallel tubes)"  annotation(Dialog(tab="General", group="Heat transfer"));
    parameter Boolean dynamicTerm=false 
      " = true, include dynamic term in momentum balance, only if not lumped_dp and not static"
                                                                                                annotation(Evaluate=true, Dialog(tab="Advanced", group="Conservation of mass, energy, momentum"));
    
    //Pressure drop and heat transfer    
    replaceable package WallFriction = 
        PressureLosses.Utilities.WallFriction.QuadraticTurbulent extends 
      Modelica_Fluid.PressureLosses.Utilities.WallFriction.PartialWallFriction 
      "Characteristic of wall friction"                                                            annotation(choicesAllMatching, Dialog(tab="General", group="Pressure drop"));
    parameter SI.Length roughness_1=2.5e-5 
      "Absolute roughness of pipe (default = smooth steel pipe)" annotation(Dialog(tab="General", group="Pressure drop"));
    parameter SI.Length roughness_2=2.5e-5 
      "Absolute roughness of pipe (default = smooth steel pipe)" annotation(Dialog(tab="General", group="Pressure drop"));
    parameter SI.DynamicViscosity eta_nominal_M1=0.01 
      "Nominal dynamic viscosity of medium 1(e.g. eta_liquidWater = 1e-3, eta_air = 1.8e-5)"
                                                                                             annotation(Dialog(tab="General", group="Pressure drop"));
    parameter SI.DynamicViscosity eta_nominal_M2=0.01 
      "Nominal dynamic viscosity of medium 1(e.g. eta_liquidWater = 1e-3, eta_air = 1.8e-5)"
                                                                                         annotation(Dialog(tab="General", group="Pressure drop"));
    parameter Boolean use_eta_nominal=false 
      "= true, if eta_nominal is used, otherwise computed from medium" annotation(Evaluate=true, Dialog(tab="General", group="Pressure drop"));
    replaceable model HeatTransfer_1 = 
        Modelica_Fluid.HeatTransfer.PipeHT_constAlpha 
        extends Modelica_Fluid.HeatTransfer.PartialPipeHeatTransfer 
      "Heat transfer model"                                                                         annotation(choicesAllMatching, Dialog(tab="General", group="Heat transfer"));
    replaceable model HeatTransfer_2 = 
        Modelica_Fluid.HeatTransfer.PipeHT_constAlpha 
        extends Modelica_Fluid.HeatTransfer.PartialPipeHeatTransfer 
      "Heat transfer model"                                                                         annotation(choicesAllMatching, Dialog(tab="General", group="Heat transfer"));
    //Display variables
    SI.HeatFlowRate Q_flow_1 "Total heat flow rate of inner pipe";
    SI.HeatFlowRate Q_flow_2 "Total heat flow rate of outer pipe";
    
    Modelica_Fluid.Components.LongPipe pipe_1(
      redeclare package Medium = Medium_1,
      n=n,
      static=static,
      lumped_dp=lumped_dp,
      kineticTerm=kineticTerm,
      crossSectionType=Modelica_Fluid.Types.CrossSectionTypes.Circular,
      length=length,
      use_wall=true,
      area_h=Modelica.Constants.pi*di_1*length*K1,
      redeclare HeatTransfer_1 heat,
      initOption=initOption_1,
      use_T_start=use_T_start_1,
      T_start=T_start_1,
      h_start=h_start_1,
      X_start=X_start_1,
      mflow_start=mflow_start_1,
      d_inner=di_1,
      d_outer=da_1,
      redeclare package WallFriction = WallFriction,
      roughness=roughness_1,
      use_eta_nominal=use_eta_nominal,
      eta_nominal=eta_nominal_M1,
      redeclare final model Wall = 
          Modelica_Fluid.Components.Wall_constProps (
          d_wall=d_wall,
          c_wall=c_wall,
          initOption=if initWall_steadyState then 3 else 2,
          T_start=T_start_wall),
        dynamicTerm=dynamicTerm) 
                               annotation (extent=[-40,-60; 20,0]);
    
    Modelica_Fluid.Components.LongPipe pipe_2(
      redeclare package Medium = Medium_2,
      n=n,
      static=static,
      lumped_dp=lumped_dp,
      kineticTerm=kineticTerm,
      crossSectionType=Modelica_Fluid.Types.CrossSectionTypes.General,
      length=length,
      redeclare HeatTransfer_2 heat,
      use_T_start=use_T_start_2,
      T_start=T_start_2,
      h_start=h_start_2,
      X_start=X_start_2,
      initOption=initOption_2,
      mflow_start=mflow_start_2,
      P_inner=Modelica.Constants.pi*(da_1 + di_2),
      A_inner=Modelica.Constants.pi/4*(di_2*di_2 - da_1*da_1),
      area_h=Modelica.Constants.pi*da_1*length*K2,
      p_a_start=p_a_start1,
      p_b_start=p_b_start2,
      redeclare package WallFriction = WallFriction,
      roughness=roughness_2,
      use_eta_nominal=use_eta_nominal,
      eta_nominal=eta_nominal_M2,
      dynamicTerm=dynamicTerm) 
                annotation (extent=[-40,88; 20,28]);
    annotation (Diagram(Line(points=[-10,36; -10,-8], style(
            color=1,
            rgbcolor={255,0,0},
            thickness=2))), Icon(
        Rectangle(extent=[-100,-26; 100,-30], style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=10,
            rgbfillColor={95,95,95},
            fillPattern=7)),
        Rectangle(extent=[-100,30; 100,26], style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=10,
            rgbfillColor={95,95,95},
            fillPattern=7)),
        Rectangle(extent=[-100,60; 100,30], style(
            color=0,
            rgbcolor={0,0,0},
            gradient=2,
            fillColor=70,
            rgbfillColor={0,63,125})),
        Rectangle(extent=[-100,-30; 100,-60], style(
            color=0,
            rgbcolor={0,0,0},
            gradient=2,
            fillColor=70,
            rgbfillColor={0,63,125})),
        Rectangle(extent=[-100,26; 100,-26], style(
            color=69,
            rgbcolor={0,128,255},
            gradient=2,
            fillColor=69,
            rgbfillColor={0,128,255})),
        Text(
          extent=[-100,-60; 100,-100],
          string="%name",
          style(color=3, rgbcolor={0,0,255}))));
    Modelica_Fluid.Interfaces.FluidPort_b port_b1(redeclare package Medium = 
          Medium_1) annotation (extent=[100,-12; 120,8]);
    Modelica_Fluid.Interfaces.FluidPort_a port_a1(redeclare package Medium = 
          Medium_1) annotation (extent=[-120,-12; -100,8]);
    Modelica_Fluid.Interfaces.FluidPort_a port_a2(redeclare package Medium = 
          Medium_2) annotation (extent=[-120,36; -100,56]);
    Modelica_Fluid.Interfaces.FluidPort_b port_b2(redeclare package Medium = 
          Medium_2) annotation (extent=[100,-56; 120,-36]);
    
  equation 
    Q_flow_1=sum(pipe_1.Qs_flow);
    Q_flow_2=sum(pipe_2.Qs_flow);
    connect(pipe_2.thermalPort, pipe_1.thermalPort);
    connect(pipe_2.port_b, port_b2) annotation (points=[20,58; 60,58; 60,-46; 110,
          -46], style(
        color=69,
        rgbcolor={0,127,255},
        thickness=2,
        gradient=2,
        fillColor=42,
        rgbfillColor={213,0,0}));
    connect(pipe_1.port_b, port_b1) annotation (points=[20,-30; 42,-30; 42,-2;
          110,-2], style(
        color=69,
        rgbcolor={0,127,255},
        thickness=2,
        gradient=2,
        fillColor=42,
        rgbfillColor={213,0,0}));
    connect(pipe_1.port_a, port_a1) annotation (points=[-40.6,-30; -75.3,-30;
            -75.3,-2; -110,-2],
                              style(
        color=69,
        rgbcolor={0,127,255},
        thickness=2,
        gradient=2,
        fillColor=42,
        rgbfillColor={213,0,0}));
    connect(pipe_2.port_a, port_a2) annotation (points=[-40.6,58; -76,58; -76,46;
            -110,46],
                    style(
        color=69,
        rgbcolor={0,127,255},
        thickness=2,
        gradient=2,
        fillColor=42,
        rgbfillColor={213,0,0}));
  end HeatExchanger;
  Sources.FixedAmbient source(
    redeclare package Medium = 
        Modelica.Media.Water.ConstantPropertyLiquidWater,
    use_T=true,
    T=Modelica.SIunits.Conversions.from_degC(20)) 
    annotation (extent=[-100,-90; -80,-70]);
  
  Components.PressureDropPipe pipeFriction(
    redeclare package Medium = 
        Modelica.Media.Water.ConstantPropertyLiquidWater,
    frictionType=Modelica_Fluid.Types.FrictionTypes.ConstantTurbulent,
    m_flow_nominal=1000,
    dp_nominal=Modelica.SIunits.Conversions.from_bar(0.3),
    flowDirection= Modelica_Fluid.Types.FlowDirection.Unidirectional) 
    annotation (extent=[-32,-40; -16,-20]);
  
  Components.StaticHead pipeHead(H_b_a=50, redeclare package Medium = 
        Modelica.Media.Water.ConstantPropertyLiquidWater,
    flowDirection= Modelica_Fluid.Types.FlowDirection.Unidirectional) 
    annotation (extent=[-48,-64; -30,-38], rotation=90);
  Components.Pump pumps(
    checkValve=true,
    redeclare package Medium = 
        Modelica.Media.Water.ConstantPropertyLiquidWater,
    N_nom=1200,
    redeclare function flowCharacteristic = 
        Modelica_Fluid.Types.PumpCharacteristics.quadraticFlow (q_nom={0,
            0.25,0.5}, head_nom={100,60,0}),
    Np_nom=4,
    M=50,
    T_start=Modelica.SIunits.Conversions.from_degC(20)) 
    annotation (extent=[-74,-88; -48,-62]);
  
  Components.Tank reservoir(
    H0=18,
    initOption=Modelica_Fluid.Types.Init.InitialValues,
    redeclare package Medium = 
        Modelica.Media.Water.ConstantPropertyLiquidWater,
    area=500,
    level_start=1.8,
    T_start=Modelica.SIunits.Conversions.from_degC(20),
    pipeArea=0.1) 
    annotation (extent=[-14,-16; 6,4]);
  
  Components.ValveLinear userValve(redeclare package Medium = 
        Modelica.Media.Water.ConstantPropertyLiquidWater, Kv=200/2e5,
    flowDirection= Modelica_Fluid.Types.FlowDirection.Unidirectional) 
    annotation (extent=[58,-38; 74,-22]);
  Sources.FixedAmbient ambient(redeclare package Medium = 
        Modelica.Media.Water.ConstantPropertyLiquidWater) 
    annotation (extent=[100,-40; 80,-20]);
  Modelica.Blocks.Sources.Step valveOpening(height=0, offset=1) 
    annotation (extent=[64,0; 84,20]);
  Modelica.Blocks.Sources.Constant PressureSetPoint(k=2e5) 
    annotation (extent=[-100,60; -80,80]);
  Modelica.Blocks.Logical.OnOffController controller(bandwidth=1000,
      pre_y_start=true) annotation (extent=[-40,60; -20,80]);
  Modelica.Blocks.Logical.TriggeredTrapezoid PumpRPMGenerator(
    rising=3,
    falling=3,
    amplitude=1200,
    offset=0.001) annotation (extent=[0,60; 20,80]);
  Sensors.RelativePressure reservoirPressure(redeclare package Medium = 
        Modelica.Media.Water.ConstantPropertyLiquidWater) 
    annotation (extent=[8,-12; 28,-32]);
  Modelica.Blocks.Continuous.FirstOrder PT1(T=50) 
    annotation (extent=[40,60; 60,80]);
  inner Components.FluidOptions fluidOptions 
    annotation (extent=[80,-100; 100,-80]);
equation 
  
  annotation (
    Diagram,
    Documentation(info="<html>
Water is pumped from a source by 4 pumps in parallel (fitted with check valves), through a pipe whose outlet is 50 m higher than the source, into a reservoir placed on a 18-m high tower. The users are represented by an equivalent valve, connected to the reservoir.
<p>
The water controller is a simple on-off controller, acting on the gauge pressure measured at the base of the tower; the output of the controller is the rotational speed of the pumps. A small but nonzero rotational speed is used to represent the standby state of the pumps, in order to avoid singularities in the flow characteristic.
<p>
Simulate for 2000 s. The pump turns on and off to keep the reservoir level around 2.5 meters, which means 20.5 meters higher than the base of the tower, corresponding to a gauge pressure of 2 bar
<p>
If using Dymola, turn off \"Equidistant time grid\" to avoid numerical errors.
</html>", revisions="<html>
<ul>
<li><i>2 Nov 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Created.</li>
</ul>
</html>"),
    experiment(StopTime=2000, NumberOfIntervals=5000),
    experimentSetupOutput(equdistant=false));
  connect(pumps.outlet, pipeHead.port_a) annotation (points=[-53.2,-70.84;
        -46,-72; -38,-72; -38,-64; -39,-64], style(color=69, rgbcolor={0,127,
          255}));
  connect(pipeFriction.port_b, reservoir.port) annotation (points=[-16,-30;
        -4,-30; -4,-17], style(color=69, rgbcolor={0,127,255}));
  connect(reservoir.port, userValve.port_a) annotation (points=[-4,-17; -4,
        -30; 58,-30],
                   style(color=69, rgbcolor={0,127,255}));
  connect(userValve.port_b, ambient.port)  annotation (points=[74,-30; 80,
        -30],
      style(color=69, rgbcolor={0,127,255}));
  connect(source.port, pumps.inlet) annotation (points=[-80,-80; -73.2,-80;
        -73.2,-77.6; -71.4,-77.6], style(color=69, rgbcolor={0,127,255}));
  connect(valveOpening.y, userValve.opening) annotation (points=[85,10; 98,
        10; 98,-12; 66,-12; 66,-22.8],
                                   style(color=74, rgbcolor={0,0,127}));
  connect(pipeHead.port_b, pipeFriction.port_a) annotation (points=[-39,-38;
        -39,-30; -32,-30],   style(color=69, rgbcolor={0,127,255}));
  connect(PressureSetPoint.y, controller.reference) annotation (points=[-79,70;
        -60,70; -60,76; -42,76], style(color=74, rgbcolor={0,0,127}));
  connect(controller.y, PumpRPMGenerator.u) 
    annotation (points=[-19,70; -2,70], style(color=5, rgbcolor={255,0,255}));
  connect(reservoir.port, reservoirPressure.port_a) annotation (points=[-4,-17;
        3,-17; 3,-22; 7,-22], style(
      color=3,
      rgbcolor={0,0,255},
      pattern=3));
  connect(reservoirPressure.p_rel, controller.u) annotation (points=[18,-13; 18,
        50; -52,50; -52,64; -42,64], style(color=74, rgbcolor={0,0,127}));
  connect(reservoirPressure.port_b, ambient.port) annotation (points=[29,-22;
        44,-22; 44,-48; 80,-48; 80,-30],
                                      style(
      color=69,
      rgbcolor={0,127,255},
      pattern=3));
  connect(PumpRPMGenerator.y, PT1.u) 
    annotation (points=[21,70; 38,70], style(color=74, rgbcolor={0,0,127}));
  connect(PT1.y, pumps.N_in) annotation (points=[61,70; 74,70; 74,30; -64.38,30;
        -64.38,-69.28], style(color=74, rgbcolor={0,0,127}));
end PumpingSystem;
