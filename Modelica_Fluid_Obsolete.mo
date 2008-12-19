import "D:/otter/Modelica/subversion/svn.Modelica.org/Modelica_Fluid/Modelica_Fluid/package.mo";

package Modelica_Fluid_Obsolete 
  package Tanks 
    model TwoTanksSimpleWater 
      import Modelica.SIunits.Conversions.*;
      import Modelica_Fluid;
      extends Modelica.Icons.Example;
      annotation (
        Diagram,
        Coordsys(grid=[1, 1], component=[20, 20]),
        experiment(StopTime=50),
        Documentation(info="<html>
<p>
Two tanks with different initial levels are connected by a 
horizontal pipe. After about 30 s, the levels of the
two tanks are identical and the system is in steady state.
The two tanks have different initial temperatures. Since the
water is flowing from tank1 to tank2, the temperature of tank2
changes until water from tank1 is flowing into tank2.
</p>
</html>"),
        experimentSetupOutput);
      Modelica_Fluid.WorkInProgress.Components.Tank Tank1(
        area=1,
        T_start=from_degC(50),
        level_start=3,
        redeclare package Medium = 
            Modelica.Media.Water.ConstantPropertyLiquidWater,
        initType=Modelica_Fluid.Types.Init.InitialValues,
        V0=0.1,
        pipeArea=0.01) 
        annotation (extent=[-70,20; -50,40]);
      Modelica_Fluid.WorkInProgress.Components.Tank Tank2(
        area=1,
        T_start=from_degC(100),
        level_start=1,
        redeclare package Medium = 
            Modelica.Media.Water.ConstantPropertyLiquidWater,
        initType=Modelica_Fluid.Types.Init.InitialValues,
        V0=0.1,
        pipeArea=0.01) 
        annotation (extent=[10,20; 30,40]);
      Modelica_Fluid.PressureLosses.PressureDropPipe shortPipe1(
        m_flow_nominal=2000,
        dp_nominal=from_bar(0.1),
        redeclare package Medium = 
            Modelica.Media.Water.ConstantPropertyLiquidWater,
        frictionType=Modelica_Fluid.Types.FrictionTypes.ConstantLaminar) 
        annotation (extent=[-29,0; -9,20]);
      
      inner Modelica_Fluid.Ambient ambient 
                                       annotation (extent=[46,64; 66,84]);
    equation 
      connect(Tank1.port, shortPipe1.port_a) 
        annotation (points=[-60,19; -60,10; -29,10],      style(color=69));
      connect(shortPipe1.port_b, Tank2.port) annotation (points=[-9,10; 20,10; 20,
            19], style(color=69, rgbcolor={0,127,255}));
    end TwoTanksSimpleWater;

    model ThreeTanksSimpleWater 
      import Modelica.SIunits.Conversions.*;
      import Modelica_Fluid;
      extends Modelica.Icons.Example;
      annotation (
        Diagram,
        Coordsys(grid=[1, 1], component=[20, 20]),
        experiment(StopTime=50),
        Documentation(info="<html>
<p>
This example demonstrates the mixing of simple (constant
property) liquid water model
between three tanks with different temperature.
The water model is incompressible.
The same tank system with a compressible water model
is provided in ThreeTanksIF97.
</p>
</html>"),
        experimentSetupOutput);
      Modelica_Fluid.WorkInProgress.Components.Tank Tank1(
        area=1,
        T_start=from_degC(50),
        level_start=3,
        redeclare package Medium = 
            Modelica.Media.Water.ConstantPropertyLiquidWater,
        initType=Modelica_Fluid.Types.Init.InitialValues,
        V0=0.1,
        pipeArea=0.01) 
        annotation (extent=[-90,20; -70,40]);
      Modelica_Fluid.WorkInProgress.Components.Tank Tank2(
        area=1,
        T_start=from_degC(100),
        level_start=1,
        redeclare package Medium = 
            Modelica.Media.Water.ConstantPropertyLiquidWater,
        initType=Modelica_Fluid.Types.Init.InitialValues,
        H0=0.5,
        V0=0.1,
        pipeArea=0.01) 
        annotation (extent=[-10, 20; 10, 40]);
      Modelica_Fluid.PressureLosses.PressureDropPipe shortPipe1(
        m_flow_nominal=2000,
        dp_nominal=from_bar(0.1),
        redeclare package Medium = 
            Modelica.Media.Water.ConstantPropertyLiquidWater,
        frictionType=Modelica_Fluid.Types.FrictionTypes.ConstantLaminar) 
        annotation (extent=[-50, -30; -30, -10]);
      Modelica_Fluid.WorkInProgress.Components.Tank Tank3(
        area=1,
        T_start=from_degC(20),
        level_start=2,
        redeclare package Medium = 
            Modelica.Media.Water.ConstantPropertyLiquidWater,
        initType=Modelica_Fluid.Types.Init.InitialValues,
        V0=0.1,
        pipeArea=0.01) 
        annotation (extent=[70,20; 90,40]);
      Modelica_Fluid.PressureLosses.PressureDropPipe shortPipe3(
        m_flow_nominal=2000,
        dp_nominal=from_bar(0.1),
        redeclare package Medium = 
            Modelica.Media.Water.ConstantPropertyLiquidWater,
        frictionType=Modelica_Fluid.Types.FrictionTypes.ConstantLaminar) 
        annotation (extent=[30, -30; 50, -10]);
      Modelica_Fluid.PressureLosses.PressureDropPipe shortPipe2(
        m_flow_nominal=1000,
        dp_nominal=from_bar(0.1),
        redeclare package Medium = 
            Modelica.Media.Water.ConstantPropertyLiquidWater,
        frictionType=Modelica_Fluid.Types.FrictionTypes.ConstantLaminar,
        medium_b(T(stateSelect=StateSelect.avoid)),
        medium_a(T(stateSelect=StateSelect.avoid))) 
        annotation (extent=[-10, -10; 10, 10], rotation=-90);
      
      inner Modelica_Fluid.Ambient ambient 
                                       annotation (extent=[49,65; 69,85]);
    equation 
      connect(Tank1.port, shortPipe1.port_a) 
        annotation (points=[-80,19; -80,-20; -50,-20],    style(color=69));
      connect(shortPipe3.port_b, Tank3.port) 
        annotation (points=[50,-20; 80,-20; 80,19],    style(color=69));
      connect(Tank2.port, shortPipe2.port_a) 
        annotation (points=[0,19; 0,10; -6.12303e-016,10],    style(color=69));
      connect(shortPipe1.port_b, shortPipe3.port_a) 
        annotation (points=[-30,-20; 30,-20],   style(color=69));
      connect(shortPipe2.port_b, shortPipe3.port_a) annotation (points=[
            6.12303e-016,-10; 0,-10; 0,-20; 30,-20],     style(color=69));
    end ThreeTanksSimpleWater;

    model ThreeTanksIF97 
      import Modelica.SIunits.Conversions.*;
      import Modelica_Fluid;
      extends Modelica.Icons.Example;
      annotation (
        Diagram,
        Coordsys(grid=[1, 1], component=[20, 20]),
        experiment(StopTime=50),
        Documentation(info="<html>
<p>
This example is the same as the \"ThreeTanksSimpleWater\"
model. The only difference is that the very detailed,
compressible medium model WaterIF97 is used.
</p>
</html>"),
        experimentSetupOutput);
      Modelica_Fluid.WorkInProgress.Components.Tank Tank1(
        area=1,
        T_start=from_degC(50),
        level_start=3,
        redeclare package Medium = Modelica.Media.Water.WaterIF97_ph,
        initType=Modelica_Fluid.Types.Init.InitialValues,
        pipeArea=0.01) 
        annotation (extent=[-90, 20; -70, 40]);
      Modelica_Fluid.WorkInProgress.Components.Tank Tank2(
        area=1,
        level_start=1,
        redeclare package Medium = Modelica.Media.Water.WaterIF97_ph,
        initType=Modelica_Fluid.Types.Init.InitialValues,
        T_start=from_degC(90),
        pipeArea=0.01) 
        annotation (extent=[-10,20; 10,40]);
      Modelica_Fluid.PressureLosses.PressureDropPipe shortPipe1(
        m_flow_nominal=2000,
        dp_nominal=from_bar(0.1),
        frictionType=Modelica_Fluid.Types.FrictionTypes.ConstantLaminar,
        redeclare package Medium = Modelica.Media.Water.WaterIF97_ph) 
        annotation (extent=[-50, -30; -30, -10]);
      Modelica_Fluid.WorkInProgress.Components.Tank Tank3(
        area=1,
        T_start=from_degC(20),
        level_start=2,
        redeclare package Medium = Modelica.Media.Water.WaterIF97_ph,
        initType=Modelica_Fluid.Types.Init.InitialValues,
        pipeArea=0.01) 
        annotation (extent=[70, 20; 90, 40]);
      Modelica_Fluid.PressureLosses.PressureDropPipe shortPipe3(
        m_flow_nominal=2000,
        dp_nominal=from_bar(0.1),
        frictionType=Modelica_Fluid.Types.FrictionTypes.ConstantLaminar,
        redeclare package Medium = Modelica.Media.Water.WaterIF97_ph) 
        annotation (extent=[30, -30; 50, -10]);
      Modelica_Fluid.PressureLosses.PressureDropPipe shortPipe2(
        m_flow_nominal=1000,
        dp_nominal=from_bar(0.1),
        frictionType=Modelica_Fluid.Types.FrictionTypes.ConstantLaminar,
        redeclare package Medium = Modelica.Media.Water.WaterIF97_ph) 
        annotation (extent=[-10,-10; 10,10],   rotation=-90);
      
      inner Modelica_Fluid.Ambient ambient 
                                       annotation (extent=[50,61; 70,81]);
    equation 
      connect(Tank1.port, shortPipe1.port_a) 
        annotation (points=[-80,19; -80,-20; -50,-20],    style(color=69));
      connect(shortPipe3.port_b, Tank3.port) 
        annotation (points=[50,-20; 80,-20; 80,19],    style(color=69));
      connect(Tank2.port, shortPipe2.port_a) 
        annotation (points=[0,19; 0,10; -6.12303e-016,10],    style(color=69));
      connect(shortPipe1.port_b, shortPipe3.port_a) 
        annotation (points=[-30,-20; 30,-20],   style(color=69));
      connect(shortPipe2.port_b, shortPipe3.port_a) annotation (points=[
            6.12303e-016,-10; 0,-10; 0,-20; 30,-20],     style(color=69));
    end ThreeTanksIF97;

    model ThreeTanksIF97SteadyState 
      import Modelica.SIunits.Conversions.*;
      import Modelica_Fluid;
      extends Modelica.Icons.Example;
      annotation (
        Diagram,
        Coordsys(grid=[1, 1], component=[20, 20]),
        experiment(StopTime=50),
        Documentation(info="<html>
<p>
This example is the same as ThreeTanks97. The only difference
is that the system starts in steady state, i.e., with constant
(mixing) temperature.
</p>
</html>"),
        experimentSetupOutput);
      Modelica_Fluid.WorkInProgress.Components.Tank Tank1(
        area=1,
        T_start=from_degC(50),
        level_start=3,
        redeclare package Medium = Modelica.Media.Water.WaterIF97_ph,
        pipeArea=0.01) 
        annotation (extent=[-90, 20; -70, 40]);
      Modelica_Fluid.WorkInProgress.Components.Tank Tank2(
        area=1,
        level_start=1,
        redeclare package Medium = Modelica.Media.Water.WaterIF97_ph,
        T_start=from_degC(90),
        pipeArea=0.01) 
        annotation (extent=[-10,20; 10,40]);
      Modelica_Fluid.PressureLosses.PressureDropPipe shortPipe1(
        m_flow_nominal=2000,
        dp_nominal=from_bar(0.1),
        frictionType=Modelica_Fluid.Types.FrictionTypes.ConstantLaminar,
        redeclare package Medium = Modelica.Media.Water.WaterIF97_ph) 
        annotation (extent=[-50, -30; -30, -10]);
      Modelica_Fluid.WorkInProgress.Components.Tank Tank3(
        area=1,
        T_start=from_degC(20),
        level_start=2,
        redeclare package Medium = Modelica.Media.Water.WaterIF97_ph,
        pipeArea=0.01) 
        annotation (extent=[70, 20; 90, 40]);
      Modelica_Fluid.PressureLosses.PressureDropPipe shortPipe3(
        m_flow_nominal=2000,
        dp_nominal=from_bar(0.1),
        frictionType=Modelica_Fluid.Types.FrictionTypes.ConstantLaminar,
        redeclare package Medium = Modelica.Media.Water.WaterIF97_ph) 
        annotation (extent=[30, -30; 50, -10]);
      Modelica_Fluid.PressureLosses.PressureDropPipe shortPipe2(
        m_flow_nominal=1000,
        dp_nominal=from_bar(0.1),
        frictionType=Modelica_Fluid.Types.FrictionTypes.ConstantLaminar,
        redeclare package Medium = Modelica.Media.Water.WaterIF97_ph) 
        annotation (extent=[-10, -10; 10, 10], rotation=-90);
      inner Modelica_Fluid.Ambient ambient 
                                       annotation (extent=[60,65; 80,85]);
    equation 
      connect(Tank1.port, shortPipe1.port_a) 
        annotation (points=[-80,19; -80,-20; -50,-20],    style(color=69));
      connect(shortPipe3.port_b, Tank3.port) 
        annotation (points=[50,-20; 80,-20; 80,19],    style(color=69));
      connect(Tank2.port, shortPipe2.port_a) 
        annotation (points=[0,19; 0,10; -6.12303e-016,10],    style(color=69));
      connect(shortPipe1.port_b, shortPipe3.port_a) 
        annotation (points=[-30,-20; 30,-20],   style(color=69));
      connect(shortPipe2.port_b, shortPipe3.port_a) annotation (points=[
            6.12303e-016,-10; 0,-10; 0,-20; 30,-20],     style(color=69));
    end ThreeTanksIF97SteadyState;

    model ThreeTanksIF97WithShortPipe2 
      import Modelica.SIunits.Conversions;
      import Modelica_Fluid;
      extends Modelica.Icons.Example;
      
      replaceable package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater 
          extends Modelica.Media.Interfaces.PartialMedium 
        "Medium in the component" 
          annotation (choicesAllMatching = true);
      // replaceable package Medium = Modelica.Media.Water.WaterIF97_ph 
      replaceable package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent;
      
      parameter Real D=0.1 "pipe diameters";
      final parameter Real A = Modelica.Constants.pi*(D/2)^2 
        "pipe area (for tank)";
      
      annotation (structurallyIncomplete,
        Diagram,
        Coordsys(grid=[1, 1], component=[20, 20]),
        experiment(StopTime=50),
        Documentation(info="<html>
<p>
Demonstrate usage of ShortPipe model with three tanks and three short pipes.
Results when using ConstantPropertyLiquidWater:
</p>
 
<img src=\"../Images/Examples/ThreeTanksResult3.png\">
 
<p>
Results when using Modelica.Media.Water.WaterIF97_ph:
</p>
 
<img src=\"../Images/Examples/ThreeTanksResult3.png\">
 
</html>"),
        experimentSetupOutput);
      Modelica_Fluid.WorkInProgress.Components.Tank Tank1(
        area=1,
        T_start=Conversions.from_degC(50),
        level_start=3,
        redeclare package Medium = Medium,
        pipeArea=A) 
        annotation (extent=[-90, 20; -70, 40]);
      Modelica_Fluid.WorkInProgress.Components.Tank Tank2(
        area=1,
        level_start=1,
        initType=Modelica_Fluid.Types.Init.InitialValues,
        T_start=Conversions.from_degC(90),
        redeclare package Medium = Medium,
        pipeArea=A) 
        annotation (extent=[-10,20; 10,40]);
      Modelica_Fluid.WorkInProgress.Components.Tank Tank3(
        area=1,
        T_start=Conversions.from_degC(20),
        level_start=2,
        redeclare package Medium = Medium,
        pipeArea=A) 
        annotation (extent=[70, 20; 90, 40]);
      Modelica_Fluid_Obsolete.ShortPipe2 shortPipe1(
        redeclare package Medium = Medium,
        length=1,
        height_ab=-1,
        T_start=Conversions.from_degC(20),
        diameter=D,
        redeclare package WallFriction = WallFriction,
        initVolume2=Modelica_Fluid.Types.Init.NoInit) 
        annotation (extent=[-50, -30; -30, -10]);
      Modelica_Fluid_Obsolete.ShortPipe2 shortPipe2(
        redeclare package Medium = Medium,
        T_start=Conversions.from_degC(20),
        length=1,
        diameter=D,
        height_ab=-1,
        redeclare package WallFriction = WallFriction,
        initVolume2=Modelica_Fluid.Types.Init.NoInit) 
        annotation (extent=[-10,-10; 10,10],   rotation=-90);
      Modelica_Fluid_Obsolete.ShortPipe2 shortPipe3(
        redeclare package Medium = Medium,
        T_start=Conversions.from_degC(20),
        length=1,
        height_ab=1,
        diameter=D,
        redeclare package WallFriction = WallFriction) 
        annotation (extent=[30, -30; 50, -10]);
      inner Modelica_Fluid.Ambient ambient annotation (extent=[-70,70; -50,90]);
    equation 
      connect(Tank1.port, shortPipe1.port_a) 
        annotation (points=[-80,19; -80,-20; -51,-20],    style(color=69));
      connect(shortPipe3.port_b, Tank3.port) 
        annotation (points=[51,-20; 80,-20; 80,19],    style(color=69));
      connect(Tank2.port, shortPipe2.port_a) 
        annotation (points=[0,19; 0,11; -6.73533e-016,11],    style(color=69));
      connect(shortPipe1.port_b, shortPipe3.port_a) 
        annotation (points=[-29,-20; 29,-20],   style(color=69));
      connect(shortPipe2.port_b, shortPipe3.port_a) annotation (points=[
            6.73533e-016,-11; 0,-11; 0,-20; 29,-20],     style(color=69));
    end ThreeTanksIF97WithShortPipe2;
  end Tanks;
  annotation (uses(Modelica(version="2.2.1"), Modelica_Fluid(version=
            "1.0 Beta 3")));
  model ShortPipe2 
    "Short pipe with two volumes, wall friction and gravity effect" 
    
    extends Modelica_Fluid.WorkInProgress.Interfaces.PartialTwoPortTransport;
    replaceable package WallFriction = 
      Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent
      extends 
      Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.PartialWallFriction
      "Characteristic of wall friction"  annotation(choicesAllMatching=true);
    
    parameter Modelica.SIunits.Length length "Length of pipe";
    parameter Modelica.SIunits.Diameter diameter 
      "Inner (hydraulic) diameter of pipe";
    parameter Modelica.SIunits.Length height_ab=0.0 
      "Height of port_b over port_a"                                   annotation(Evaluate=true);
    parameter Modelica.SIunits.Length roughness(min=0)=2.5e-5 
      "Absolute roughness of pipe (default = smooth steel pipe)" 
        annotation(Dialog(enable=WallFriction.use_roughness));
    parameter Boolean use_nominal = false 
      "= true, if eta_nominal and d_nominal are used, otherwise computed from medium"
        annotation(Evaluate=true);
    parameter Modelica.SIunits.DynamicViscosity eta_nominal=0.01 
      "Nominal dynamic viscosity (for wall friction computation)" annotation(Dialog(enable=use_nominal));
    parameter Modelica.SIunits.Density d_nominal=0.01 
      "Nominal density (for wall friction computation)" annotation(Dialog(enable=use_nominal));
  /*
  parameter Modelica_Fluid.Types.FlowDirectionWithGlobalDefault.Temp 
    flowDirection=
            Modelica_Fluid.Types.FlowDirectionWithGlobalDefault.UseGlobalFluidOption 
    "Unidirectional (port_a -> port_b) or bidirectional flow component" 
     annotation(Dialog(tab="Advanced"));
*/
    parameter Modelica.SIunits.AbsolutePressure dp_small=1 
      "Turbulent flow for wall friction if |dp| >= dp_small" 
      annotation(Dialog(tab="Advanced", enable=WallFriction.use_dp_small));
    final parameter Modelica.SIunits.Volume V=Modelica.Constants.pi*(diameter/2)
        ^2*length;
    
    parameter Modelica_Fluid.Types.Init.Temp initVolume1=
               Modelica_Fluid.Types.Init.NoInit 
      "Initialization option for volume 1" 
      annotation(Dialog(tab = "Initialization"));
    
    parameter Modelica_Fluid.Types.Init.Temp initVolume2=
              Modelica_Fluid.Types.Init.NoInit 
      "Initialization option for volume 2" 
      annotation(Dialog(tab = "Initialization"));
    
  /*
  parameter Medium.AbsolutePressure p_start = Medium.p_default 
    "Start value of pressure" 
    annotation(Dialog(tab = "Initialization"));
  parameter Boolean use_T_start = true "Use T_start if true, otherwise h_start"
    annotation(Dialog(tab = "Initialization"), Evaluate=true);
  parameter Medium.Temperature T_start=
    if use_T_start then Medium.T_default else Medium.temperature_phX(p_start,h_start,X_start) 
    "Start value of temperature" 
    annotation(Dialog(tab = "Initialization", enable = use_T_start));
  parameter Medium.SpecificEnthalpy h_start=
    if use_T_start then Medium.specificEnthalpy_pTX(p_start, T_start, X_start) else Medium.h_default 
    "Start value of specific enthalpy" 
    annotation(Dialog(tab = "Initialization", enable = not use_T_start));
  parameter Medium.MassFraction X_start[Medium.nX] = Medium.X_default 
    "Start value of mass fractions m_i/m" 
    annotation (Dialog(tab="Initialization", enable=Medium.nXi > 0));
*/
    
    annotation (defaultComponentName="pipe",Icon(
        Rectangle(extent=[-100,60; 100,-60],   style(
            color=0,
            gradient=2,
            fillColor=8)),
        Rectangle(extent=[-100,34; 100,-36],   style(
            color=69,
            gradient=2,
            fillColor=69)),
        Text(
          extent=[-150,-60; 150,-110],
          string="%name",
          style(gradient=2, fillColor=69)),
        Ellipse(extent=[-90,15; -60,-15], style(
              color=0,
              rgbcolor={0,0,0},
              fillColor=0,
              rgbfillColor={0,0,0})),
        Ellipse(extent=[60,15; 90,-15],   style(
              color=0,
              rgbcolor={0,0,0},
              fillColor=0,
              rgbfillColor={0,0,0}))),       Documentation(info="<html>
<p>
Simple pipe model consisting of two volumes, 
wall friction (with different friction correlations)
and gravity effect. This model is mostly used to demonstrate how
to build up more detailed models from the basic components.
</p>
</html>"),
      Diagram,
      Coordsys(grid=[1,1], scale=0));
    Modelica_Fluid.PressureLosses.WallFrictionAndGravity frictionAndGravity(
      redeclare package Medium = Medium, 
      flowDirection=flowDirection, 
      redeclare package WallFriction = WallFriction, 
      diameter=diameter, 
      roughness=roughness, 
      use_nominal=use_nominal, 
      eta_nominal=eta_nominal, 
      d_nominal=d_nominal, 
      from_dp=true, 
      dp_small=dp_small, 
      show_Re=false, 
      length=length, 
      height_ab=height_ab) 
                         annotation (extent=[-10,-10; 10,10]);
    Modelica_Fluid.Pipes.BaseClasses.PortVolume volume1(
      redeclare package Medium = Medium,
      p_start=p_a_start,
      use_T_start=use_T_start,
      T_start=T_start,
      h_start=h_start,
      X_start=X_start,
      V=V/2,
      initType=initVolume1) 
      annotation (extent=[-70,-10; -50,10]);
    Modelica_Fluid.Pipes.BaseClasses.PortVolume volume2(
      redeclare package Medium = Medium,
      p_start=p_b_start,
      use_T_start=use_T_start,
      T_start=T_start,
      h_start=h_start,
      X_start=X_start,
      V=V/2,
      initType=initVolume2) 
      annotation (extent=[50,-10; 70,10]);
  equation 
    connect(volume1.port, port_a) 
      annotation (points=[-60,0; -110,0], style(color=69, rgbcolor={0,127,255}));
    connect(volume1.port, frictionAndGravity.port_a) 
      annotation (points=[-60,0; -10,0], style(color=69, rgbcolor={0,127,255}));
    connect(frictionAndGravity.port_b, volume2.port) 
      annotation (points=[10,0; 60,0], style(color=69, rgbcolor={0,127,255}));
    connect(volume2.port, port_b) 
      annotation (points=[60,0; 110,0], style(color=69, rgbcolor={0,127,255}));
  end ShortPipe2;
end Modelica_Fluid_Obsolete;
