package Subsystems 
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
    parameter Modelica_Fluid.Types.Init.Temp initType_1 "Initialization option"
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
    parameter Modelica_Fluid.Types.Init.Temp initType_2 "Initialization option"
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
        BaseClasses.PressureLosses.WallFriction.QuadraticTurbulent extends 
      Modelica_Fluid.BaseClasses.PressureLosses.WallFriction.PartialWallFriction
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
        Modelica_Fluid.BaseClasses.Pipes.HeatTransfer.PipeHT_constAlpha 
        extends 
      Modelica_Fluid.BaseClasses.Pipes.HeatTransfer.PartialPipeHeatTransfer 
      "Heat transfer model"                                                                         annotation(choicesAllMatching, Dialog(tab="General", group="Heat transfer"));
    replaceable model HeatTransfer_2 = 
        Modelica_Fluid.BaseClasses.Pipes.HeatTransfer.PipeHT_constAlpha 
        extends 
      Modelica_Fluid.BaseClasses.Pipes.HeatTransfer.PartialPipeHeatTransfer 
      "Heat transfer model"                                                                         annotation(choicesAllMatching, Dialog(tab="General", group="Heat transfer"));
    //Display variables
    SI.HeatFlowRate Q_flow_1 "Total heat flow rate of inner pipe";
    SI.HeatFlowRate Q_flow_2 "Total heat flow rate of outer pipe";
    
    Modelica_Fluid.Components.Pipes.DistributedPipeFV pipe_1(
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
      initType=initType_1,
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
          Modelica_Fluid.Components.Thermal.WallConstProps (
          d_wall=d_wall,
          c_wall=c_wall,
          initOption=if initWall_steadyState then 3 else 2,
          T_start=T_start_wall),
        dynamicTerm=dynamicTerm) 
                               annotation (extent=[-40,-60; 20,0]);
    
    Modelica_Fluid.Components.Pipes.DistributedPipeFV pipe_2(
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
      initType=initType_2,
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
          -75.3,-2; -110,-2], style(
        color=69,
        rgbcolor={0,127,255},
        thickness=2,
        gradient=2,
        fillColor=42,
        rgbfillColor={213,0,0}));
    connect(pipe_2.port_a, port_a2) annotation (points=[-40.6,58; -76,58; -76,
          46; -110,46],
                    style(
        color=69,
        rgbcolor={0,127,255},
        thickness=2,
        gradient=2,
        fillColor=42,
        rgbfillColor={213,0,0}));
  end HeatExchanger;
end Subsystems;
