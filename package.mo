package Modelica_Fluid "Fluid package that should be included into package Modelica"
annotation (
  version="0.794", 
  versionDate="2004-05-31", 
  preferedView="info", 
  Settings(NewStateSelection=true), 
  Documentation(info="<html>
<p>
This library provides basic components and property models
to model <b>1-dimensional thermo-fluid flow</b> systems of
a <b>single substance</b> or of a <b>mixture of substances</b> 
with optional <b>multiple phases</b>. The goal is to include 
this library in the Modelica standard library. The Modelica_Fluid
library uses the media models from the Modelica_Media library
</p>

<p>
The Modelica_Fluid library is still far away from a first release.
It is currently a beta release and components may be changed without
providing automatic conversion to a new version.
</p>

<p><b>Copyright &copy; 2002-2003, Modelica Association.</b></p>
<p><i>
This Modelica package is <b>free</b> software; it can be redistributed and/or modified
under the terms of the <b>Modelica license</b>, see the license conditions
and the accompanying <b>disclaimer</b> in the documentation of package
Modelica in file \"Modelica/package.mo\".
</i></p>
</html>"));


extends Modelica.Icons.Library;


package Tutorial "Tutorial" 
  annotation (DocumentationClass=true, Documentation(info="<HTML>
<h3><font color=\"#008000\" size=5>Tutorial of package Modelica_Fluid</font></h3>
<p>
Library <b>Modelica_Fluid</b> is a <b>free</b> Modelica package providing
a standardized interface to components describing
1-dimensional fluid flow in networks of pipes. A unique feature is that the
component equations and the media models are decoupled.
All components are implemented such that they can be used for
media from the Modelica_Media library. This means especially that an
incompressible or compressible medium, a single or a multiple
substance medium with one or more phases might be used for 
every component model in the Modelica_Fluid library.
</p>

<p>
This tutorial is just a start and will be improved
considerably.
</p>
</HTML>"));
  
  model ReleaseNotes "Release notes" 
    annotation (Documentation(info="<HTML>
<h3><font color=\"#008000\" size=5>Release notes</font></h3>

<h3><font color=\"#008000\">Version 0.794, 2004-05-31</font></h3>
<ul>
<li> Sensors.mo, Examples/DrumBoiler.mo: extend sensors with user choice
     for measurement unit.</li>

<li> Components.mo, Types.mo: moved components and types to 
     package Examples.</li>

<li> Moved Examples from file Modelica_Media/package.mo to own 
     Modelica_Media/Examples subdirectory and created separate 
     file per sub-package. This shall simplify the maintenance of
     examples by different authors</li>

<li> Moved Interfaces from file Modelica_Media/package.mo to 
     Modelica_Media/Interfaces.mo</li>

</ul>

<h3><font color=\"#008000\">Version 0.793, 2004-05-18</font></h3>
<ul>
<li> Removed \"semiLinear\" function since available as
     Modelica 2.1 built-in operator in Dymola.</li>

<li> Minor bug in \"Components.ShortPipe\" corrected.</li>

<li> Bug in \"Components.Orifice\" corrected
     (dp was previously calculated in
      Interfaces.PartialTwoPortTransport,
      but this was removed and not updated in Orifice).</li>
</ul>

<h3><font color=\"#008000\">Version 0.792, 2003-11-07</font></h3>
<p>
This is the first consolidated version made up from
several changes for Modelica'2003. 
Modelica_Fluid is still quite far away
from a library that could be included in the Modelica
standard library.
</p>


<h3><font color=\"#008000\">Previous Releases</font></h3>
<ul>
<li><i>Nov. 6, 2002</i><br>
       by Hilding Elmqvist: first version of the basic design.</li>
<li><i>Nov. 11, 2002</i><br>
       by Hilding Elmqvist, Martin Otter: improved version.</li>
<li><i>Nov. 20-21, 2002</i><br>
       by Hilding Elmqvist, Mike Tiller, Allan Watson, John Batteh, Chuck Newman,
       Jonas Eborn: Improved at the 32nd Modelica Design Meeting.</li>
<li><i>Nov. 30, 2002</i><br>
       by Martin Otter: Improved the design from the design meeting:
       Adapted to Modelica standard library 1.5,
       added \"choicesAllMatching=true\" annotation,
       added short documentation to \"Interfaces\",
       added packages \"Examples\" and \"Media\" (previously called \"Properties\")
       from previous versions and adapted them to the updated
       \"Interfaces\" package.</li>
<li><i>Dec., 2002</i><br>
       by Hubertus Tummescheit:
       Improved version of the high precision water model
       (Copy from ThermoFluid library, code reorganization,
       enhanced documentation, additional functions).</li>
<li><i>Feb., 2003</i><br>
       by Martin Otter: Included several elementary components and
       a model for moisted air. Some elementary components, such as
       FixedAmbient, are adapted versions from the SimpleFlow fluid library
       of Anton Haumer.</li>
<li><i>Aug., 2003</i><br>
       by Martin Otter: Improved documentation, PortVicinity (now called semiLinear)
       manually expanded, two different volume types,
       replaced number of massFractions from n to n-1 in order
       that usage of model for single substances is easier
       and in order that no special cases have to be treated
       in the equations (previously the massFraction equations had to
       be removed for single substance flow; now they are removed
       automatically, since the dimensions are zero, and not one
       as previously), included asserts to check the validity of
       the medium models, included the dynamic viscosity in the
       medium models, adapted the examples and medium models to the
       changes in Interfaces, improved menus according to the new
       features in Dymola 5.1. Added \"Components.ShortPipe\" that
       contains a detailed model of the frictional losses in pipes
       over a very wide range.</li>
<li><i>Sept., 2003</i><br>
       by Martin Otter: Changes according to the decisions of the
       Modelica design meeting in Dearborn, Sept. 2-4, 2003.
       Fluid library splitt in to two packages: Modelica_Media
       that contains the media models and Modelica_Fluid that
       contains fluid flow components. Modelica_Media is
       independent of Modelica_Fluid and my be used also from
       other packages that may have a different design as
       Modelica_Fluid.
<li><i>Oct., 2003</i><br>
       by Martin Otter: Adapted to latest design of the Modelica_Media
       library.<br>
       by Ruediger Franke: Included sensor components and 
       Modelica_Fluid.Examples.DrumBoiler example.</li>
</ul>

</HTML>
"));
  equation 
    
  end ReleaseNotes;
end Tutorial;


replaceable package PackageMedium = Modelica_Media.Interfaces.PartialMedium 
  "To allow change of default medium for all components" annotation (
    choicesAllMatching=true);

<<<<<<< package.mo

package Examples "Example models to show how to use the Fluid package" 
  import Modelica.Icons;
  
  import Modelica.Constants;
  
  import SI = Modelica.SIunits;
  
  extends Icons.Library;
  
  package TestComponents 
    "Test components (this package will be removed for the final version)" 
    
    extends Icons.Library;
    
    model TestShortPipe "Test ShortPipe componet" 
      import Modelica.SIunits.Conversions.*;
      
      extends Modelica.Icons.Example;
      annotation (
        Diagram,
        experiment(StopTime=1),
        Coordsys(grid=[1, 1], component=[20, 20]));
      Components.ShortPipe ShortPipe3(
        from_dp=false,
        dp_nominal=from_bar(0.1),
        frictionType=Modelica_Fluid.Types.FrictionTypes.DetailedFriction,
        diameter=20,
        roughness=2e-5,
        redeclare package Medium = Modelica_Media.Water.IF97_ph) 
        annotation (extent=[-2, 0; 18, 20]);
      Sources.MassFlowSource MassFlowSource3(T_ambient=from_degC(30),
          redeclare package Medium = Modelica_Media.Water.IF97_ph) 
        annotation (extent=[-41, 0; -21, 20]);
      Modelica_Fluid.Sources.FixedAmbient ambient3(T_ambient=from_degC(15),
          redeclare package Medium = Modelica_Media.Water.IF97_ph) 
        annotation (extent=[59, 0; 39, 20]);
      Modelica.Blocks.Sources.Ramp ramp1(
        height={6},
        duration={3},
        offset={-3}) annotation (extent=[-80, 0; -60, 20]);
    equation 
      connect(MassFlowSource3.port, ShortPipe3.port_a) 
        annotation (points=[-20, 10; -3, 10], style(color=69));
      connect(ShortPipe3.port_b, ambient3.port) 
        annotation (points=[19, 10; 38, 10], style(color=69));
      connect(ramp1.outPort, MassFlowSource3.m_dot) 
        annotation (points=[-59, 10; -43, 10], style(color=3));
    end TestShortPipe;
    
    model TestLongPipe "Test ShortPipe componet" 
      import Modelica.SIunits.Conversions.*;
      
      extends Modelica.Icons.Example;
      Modelica_Fluid.Sources.FixedAmbient ambient(T_ambient=from_degC(15),
          redeclare package Medium = Modelica_Media.Air.SimpleAir) 
        annotation (extent=[60, 0; 40, 20]);
      annotation (Diagram, experiment(StopTime=4));
      Components.LongPipe LongPipe(
        frictionType=Modelica_Fluid.Types.FrictionTypes.ConstantLaminar,
        redeclare package Medium = Modelica_Media.Air.SimpleAir,
        dp_nominal=from_bar(0.01),
        diameter=0.05,
        length=1) annotation (extent=[0, 0; 20, 20]);
      
      Sources.MassFlowSource MassFlowSource1(T_ambient=from_degC(30),
          redeclare package Medium = Modelica_Media.Air.SimpleAir) 
        annotation (extent=[-40, 0; -20, 20]);
      Modelica.Blocks.Sources.Ramp ramp(
        height={6},
        duration={3},
        offset={-3}) annotation (extent=[-80, 0; -60, 20]);
    equation 
      connect(LongPipe.port_b, ambient.port) 
        annotation (points=[21, 10; 39, 10], style(color=69));
      connect(MassFlowSource1.port, LongPipe.port_a) 
        annotation (points=[-19, 10; -1, 10], style(color=69));
      connect(ramp.outPort, MassFlowSource1.m_dot) 
        annotation (points=[-59, 10; -42, 10], style(color=3));
    end TestLongPipe;
    
    model TestCheckStateSelection 
      "Check whether for the choosen medium the expected states are selected" 
      
      import Modelica.SIunits.Conversions.*;
      
      Modelica_Fluid.Interfaces.JunctionVolume junctionVolume1(
        V=1.e-4,
        T_start=from_degC(50.0),
        redeclare package Medium = Modelica_Media.Water.SimpleLiquidWater) 
        annotation (extent=[-10, 90; 10, 70], rotation=180);
      extends Modelica.Icons.Example;
      Modelica_Fluid.Sources.FixedAmbient fixedAmbient1a(
        p_ambient=from_bar(1.0),
        redeclare package Medium = Modelica_Media.Water.SimpleLiquidWater,
        T_ambient=from_degC(10)) annotation (extent=[-80, 70; -60, 90]);
      Modelica_Fluid.Components.ShortPipe shortPipe1a(
        m_dot_nominal=10,
        frictionType=Modelica_Fluid.Types.FrictionTypes.ConstantLaminar,
      redeclare package Medium = Modelica_Media.Water.SimpleLiquidWater) 
        annotation (extent=[-40, 70; -20, 90]);
      annotation (
        Diagram(
          Text(
            extent=[-100, 102; -52, 92],
            string="SimpleLiquidWater",
            style(color=58)),
          Text(
            extent=[-113, 55; -54, 47],
            style(color=58),
            string="SimpleAir"),
          Text(
            extent=[-113, 5; -54, -3],
            style(color=58),
            string="DetailedAir")),
        Documentation(info="<html>
<p>
This model is used to select a medium model and then translate.
The log should give information about the selected states.
They should be the ones choosen as StateSelect.prefer in the
corresponding medium model.
</p>
</html>
"),     Coordsys(grid=[1, 1], component=[20, 20]));
      
      Modelica_Fluid.Components.ShortPipe shortPipe1b(
        m_dot_nominal=10,
        redeclare package Medium = Modelica_Media.Water.SimpleLiquidWater,
        frictionType=Modelica_Fluid.Types.FrictionTypes.ConstantLaminar) 
        annotation (extent=[20, 70; 40, 90]);
      Modelica_Fluid.Sources.FixedAmbient fixedAmbient1b(
        p_ambient=from_bar(0.95),
        T_ambient=from_degC(30),
        redeclare package Medium = Modelica_Media.Water.SimpleLiquidWater) 
        annotation (extent=[80, 70; 60, 90]);
      
      Modelica_Fluid.Interfaces.JunctionVolume junctionVolume2(
        redeclare package Medium = Modelica_Media.Air.SimpleAir,
        p_start=from_bar(0.96),
        T_start=from_degC(50.0),
        V=0.1) annotation (extent=[-11, 40; 9, 20], rotation=180);
      
      Modelica_Fluid.Sources.FixedAmbient fixedAmbient2a(
        p_ambient=from_bar(1.0),
        T_ambient=from_degC(10),
        redeclare package Medium = Modelica_Media.Air.SimpleAir) 
        annotation (extent=[-80, 20; -60, 40]);
      Modelica_Fluid.Components.ShortPipe shortPipe2a(
        redeclare package Medium = Modelica_Media.Air.SimpleAir,
        dp_nominal=from_bar(0.01),
        m_dot_nominal=1,
        frictionType=Modelica_Fluid.Types.FrictionTypes.ConstantLaminar) 
        annotation (extent=[-40, 20; -20, 40]);
      
      Modelica_Fluid.Components.ShortPipe shortPipe2b(
        redeclare package Medium = Modelica_Media.Air.SimpleAir,
        dp_nominal=from_bar(0.01),
        m_dot_nominal=1,
        frictionType=Modelica_Fluid.Types.FrictionTypes.ConstantLaminar) 
        annotation (extent=[20, 20; 40, 40]);
      
      Modelica_Fluid.Sources.FixedAmbient fixedAmbient2b(
        p_ambient=from_bar(0.95),
        redeclare package Medium = Modelica_Media.Air.SimpleAir,
        T_ambient=from_degC(10)) annotation (extent=[80, 20; 60, 40]);
      
      Modelica_Fluid.Interfaces.JunctionVolume junctionVolume3(
        p_start=from_bar(0.96),
        T_start=from_degC(50.0),
        V=0.1,
        redeclare package Medium = Modelica_Media.Air.DetailedAir) 
        annotation (extent=[-11, -10; 9, -30], rotation=180);
      Modelica_Fluid.Sources.FixedAmbient fixedAmbient3a(
        p_ambient=from_bar(1.0),
        T_ambient=from_degC(10),
        redeclare package Medium = Modelica_Media.Air.DetailedAir) 
        annotation (extent=[-80, -30; -60, -10]);
      Modelica_Fluid.Components.ShortPipe shortPipe3a(
        dp_nominal=from_bar(0.01),
        m_dot_nominal=1,
        redeclare package Medium = Modelica_Media.Air.DetailedAir,
        frictionType=Modelica_Fluid.Types.FrictionTypes.ConstantLaminar) 
        annotation (extent=[-40, -30; -20, -10]);
      Modelica_Fluid.Components.ShortPipe shortPipe3b(
        dp_nominal=from_bar(0.01),
        m_dot_nominal=1,
        redeclare package Medium = Modelica_Media.Air.DetailedAir,
        frictionType=Modelica_Fluid.Types.FrictionTypes.ConstantLaminar) 
        annotation (extent=[20, -30; 40, -10]);
      Modelica_Fluid.Sources.FixedAmbient fixedAmbient3b(
        p_ambient=from_bar(0.95),
        T_ambient=from_degC(10),
        redeclare package Medium = Modelica_Media.Air.DetailedAir) 
        annotation (extent=[80, -30; 60, -10]);
    equation 
      connect(fixedAmbient1a.port, shortPipe1a.port_a) 
        annotation (points=[-59, 80; -41, 80], style(color=69));
      connect(shortPipe1b.port_b, fixedAmbient1b.port) 
        annotation (points=[41, 80; 59, 80], style(color=69));
      connect(shortPipe1b.port_a, junctionVolume1.port) 
        annotation (points=[19, 80; 0, 80], style(color=69));
      connect(shortPipe1a.port_b, junctionVolume1.port) 
        annotation (points=[-19, 80; 0, 80], style(color=69));
      connect(fixedAmbient2a.port, shortPipe2a.port_a) 
        annotation (points=[-59, 30; -41, 30], style(color=69));
      connect(shortPipe2b.port_b, fixedAmbient2b.port) 
        annotation (points=[41, 30; 59, 30], style(color=69));
      connect(shortPipe2b.port_a, junctionVolume2.port) 
        annotation (points=[19, 30; -1, 30], style(color=69));
      connect(shortPipe2a.port_b, junctionVolume2.port) 
        annotation (points=[-19, 30; -1, 30], style(color=69));
      connect(fixedAmbient3a.port, shortPipe3a.port_a) 
        annotation (points=[-59, -20; -41, -20], style(color=69));
      connect(shortPipe3b.port_b, fixedAmbient3b.port) 
        annotation (points=[41, -20; 59, -20], style(color=69));
      connect(shortPipe3b.port_a, junctionVolume3.port) 
        annotation (points=[19, -20; -1, -20], style(color=69));
      connect(shortPipe3a.port_b, junctionVolume3.port) 
        annotation (points=[-19, -20; -1, -20], style(color=69));
    end TestCheckStateSelection;
    
    model TwoVolumesAir 
      "This example shows the difference of a JunctionVolume and of a FixedComponentVolume"
      
      
      import SI = Modelica.SIunits;
      parameter SI.Volume V=0.05 "Size of volume";
      
      Interfaces.JunctionVolume junctionVolume1(
        T_start=from_degC(50.0),
        V=V/2,
        redeclare package Medium = Modelica_Media.Air.DetailedAir,
        initType=Modelica_Media.Interfaces.PartialMedium.Choices.Init.
            InitialStates) annotation (extent=[-30, 40; -10, 20], rotation=180);
      
      import Modelica.SIunits.Conversions.*;
      
      extends Modelica.Icons.Example;
      Sources.FixedAmbient fixedAmbient1(
        p_ambient=from_bar(1.0),
        T_ambient=from_degC(10),
        redeclare package Medium = Modelica_Media.Air.DetailedAir) 
        annotation (extent=[-100, 20; -80, 40]);
      Components.ShortPipe shortPipe1(
        m_dot_nominal=10,
        frictionType=Modelica_Fluid.Types.FrictionTypes.ConstantLaminar,
        redeclare package Medium = Modelica_Media.Air.DetailedAir) 
        annotation (extent=[-60, 20; -40, 40]);
      annotation (Diagram, Documentation(info="<html>
<p>
This example shows two ways of merging flows. In the upper part
a volume in the connection port is modeled. At the initial time the
fluid in this volume has a higher temperature as the fluid
flowing form the ambient. Therefore, im shortPipe the temperature
will need some time to arrive at the mixing temperature.
</p>
<p>
In the lower part the volume in the connection port is neglected
and ideal mixing takes place. Since there is a steady state mass
flow and no (existing) delay mechanism is modeled, the mixing
is instantaneous. The temperature is the same as in the upper part.
</p>
</html>
"));
      Components.ShortPipe shortPipe3(
        m_dot_nominal=10,
        frictionType=Modelica_Fluid.Types.FrictionTypes.ConstantLaminar,
        redeclare package Medium = Modelica_Media.Air.DetailedAir) 
        annotation (extent=[40, 20; 60, 40]);
      Sources.FixedAmbient fixedAmbient3(
        T_ambient=from_degC(20),
        p_ambient=from_bar(0.95),
        redeclare package Medium = Modelica_Media.Air.DetailedAir) 
        annotation (extent=[100, 20; 80, 40]);
      Sources.FixedAmbient fixedAmbient2(
        p_ambient=from_bar(1.0),
        T_ambient=from_degC(20),
        redeclare package Medium = Modelica_Media.Air.DetailedAir) 
        annotation (extent=[-100, 60; -80, 80]);
      Components.ShortPipe shortPipe2(
        m_dot_nominal=10,
        frictionType=Modelica_Fluid.Types.FrictionTypes.ConstantLaminar,
        redeclare package Medium = Modelica_Media.Air.DetailedAir) 
        annotation (extent=[-60, 60; -40, 80]);
      Interfaces.JunctionVolume junctionVolume2(
        T_start=from_degC(50.0),
        V=V/2,
        redeclare package Medium = Modelica_Media.Air.DetailedAir) 
        annotation (extent=[0, 40; 20, 20], rotation=180);
      Interfaces.JunctionVolume junctionVolume3(
        T_start=from_degC(50.0),
        V=V,
        redeclare package Medium = Modelica_Media.Air.DetailedAir) 
        annotation (extent=[-30, -50; -10, -70], rotation=180);
      Sources.FixedAmbient fixedAmbient4(
        p_ambient=from_bar(1.0),
        T_ambient=from_degC(10),
        redeclare package Medium = Modelica_Media.Air.DetailedAir) 
        annotation (extent=[-100, -70; -80, -50]);
      Components.ShortPipe shortPipe4(
        m_dot_nominal=10,
        frictionType=Modelica_Fluid.Types.FrictionTypes.ConstantLaminar,
        redeclare package Medium = Modelica_Media.Air.DetailedAir) 
        annotation (extent=[-60, -70; -40, -50]);
      Components.ShortPipe shortPipe5(
        m_dot_nominal=10,
        frictionType=Modelica_Fluid.Types.FrictionTypes.ConstantLaminar,
        redeclare package Medium = Modelica_Media.Air.DetailedAir) 
        annotation (extent=[0, -70; 20, -50]);
      Sources.FixedAmbient fixedAmbient5(
        T_ambient=from_degC(20),
        p_ambient=from_bar(0.95),
        redeclare package Medium = Modelica_Media.Air.DetailedAir) 
        annotation (extent=[60, -70; 40, -50]);
      Sources.FixedAmbient fixedAmbient6(
        p_ambient=from_bar(1.0),
        T_ambient=from_degC(20),
        redeclare package Medium = Modelica_Media.Air.DetailedAir) 
        annotation (extent=[-100, -30; -80, -10]);
      Components.ShortPipe shortPipe6(
        m_dot_nominal=10,
        frictionType=Modelica_Fluid.Types.FrictionTypes.ConstantLaminar,
        redeclare package Medium = Modelica_Media.Air.DetailedAir) 
        annotation (extent=[-60, -30; -40, -10]);
    equation 
      connect(fixedAmbient1.port, shortPipe1.port_a) 
        annotation (points=[-79, 30; -61, 30], style(color=69));
      connect(shortPipe3.port_b, fixedAmbient3.port) 
        annotation (points=[61, 30; 79, 30], style(color=69));
      connect(fixedAmbient2.port, shortPipe2.port_a) 
        annotation (points=[-79, 70; -61, 70], style(color=69));
      connect(shortPipe2.port_b, junctionVolume1.port) 
        annotation (points=[-39, 70; -20, 70; -20, 30], style(color=69));
      connect(shortPipe1.port_b, junctionVolume1.port) 
        annotation (points=[-39, 30; -20, 30], style(color=69));
      connect(junctionVolume2.port, shortPipe3.port_a) 
        annotation (points=[10, 30; 39, 30], style(color=69));
      connect(junctionVolume1.port, junctionVolume2.port) 
        annotation (points=[-20, 30; 10, 30], style(color=69));
      connect(fixedAmbient4.port, shortPipe4.port_a) 
        annotation (points=[-79, -60; -61, -60], style(color=69));
      connect(shortPipe5.port_b, fixedAmbient5.port) 
        annotation (points=[21, -60; 39, -60], style(color=69));
      connect(shortPipe5.port_a, junctionVolume3.port) 
        annotation (points=[-1, -60; -20, -60], style(color=69));
      connect(fixedAmbient6.port, shortPipe6.port_a) 
        annotation (points=[-79, -20; -61, -20], style(color=69));
      connect(shortPipe6.port_b, junctionVolume3.port) 
        annotation (points=[-39, -20; -20, -20; -20, -60], style(color=69));
      connect(shortPipe4.port_b, junctionVolume3.port) 
        annotation (points=[-39, -60; -20, -60], style(color=69));
    end TwoVolumesAir;
    
    model TwoVolumesDetailedWater 
      "This example shows the difference of a JunctionVolume and of a FixedComponentVolume"
      
      
      import SI = Modelica.SIunits;
      parameter SI.Volume V=1.e-4 "Size of volume";
      
      Interfaces.JunctionVolume junctionVolume1(
        T_start=from_degC(50.0),
        V=V/2,
      redeclare package Medium = Modelica_Media.Water.IF97_ph) 
        annotation (extent=[-30, 40; -10, 20], rotation=180);
      
      import Modelica.SIunits.Conversions.*;
      
      extends Modelica.Icons.Example;
      Sources.FixedAmbient fixedAmbient1(
        p_ambient=from_bar(1.0),
        T_ambient=from_degC(10),
      redeclare package Medium = Modelica_Media.Water.IF97_ph) 
        annotation (extent=[-100, 20; -80, 40]);
      Components.ShortPipe shortPipe1(
        m_dot_nominal=10,
        frictionType=Modelica_Fluid.Types.FrictionTypes.ConstantLaminar,
      redeclare package Medium = Modelica_Media.Water.IF97_ph) 
        annotation (extent=[-60, 20; -40, 40]);
      annotation (Diagram, Documentation(info="<html>
<p>
This example shows two ways of merging flows. In the upper part
a volume in the connection port is modeled. At the initial time the
fluid in this volume has a higher temperature as the fluid
flowing form the ambient. Therefore, im shortPipe the temperature
will need some time to arrive at the mixing temperature.
</p>
<p>
In the lower part the volume in the connection port is neglected
and ideal mixing takes place. Since there is a steady state mass
flow and no (existing) delay mechanism is modeled, the mixing
is instantaneous. The temperature is the same as in the upper part.
</p>
</html>
"));
      Components.ShortPipe shortPipe3(
        m_dot_nominal=10,
        frictionType=Modelica_Fluid.Types.FrictionTypes.ConstantLaminar,
      redeclare package Medium = Modelica_Media.Water.IF97_ph) 
        annotation (extent=[40, 20; 60, 40]);
      Sources.FixedAmbient fixedAmbient3(
        T_ambient=from_degC(20),
        p_ambient=from_bar(0.95),
      redeclare package Medium = Modelica_Media.Water.IF97_ph) 
        annotation (extent=[100, 20; 80, 40]);
      Sources.FixedAmbient fixedAmbient2(
        p_ambient=from_bar(1.0),
        T_ambient=from_degC(20),
      redeclare package Medium = Modelica_Media.Water.IF97_ph) 
        annotation (extent=[-100, 60; -80, 80]);
      Components.ShortPipe shortPipe2(
        m_dot_nominal=10,
        frictionType=Modelica_Fluid.Types.FrictionTypes.ConstantLaminar,
      redeclare package Medium = Modelica_Media.Water.IF97_ph) 
        annotation (extent=[-60, 60; -40, 80]);
      Interfaces.JunctionVolume junctionVolume2(
        T_start=from_degC(50.0),
        V=V/2,
        initType=Modelica_Media.Interfaces.PartialMedium.Choices.Init.
            InitialStates,
      redeclare package Medium = Modelica_Media.Water.IF97_ph) 
                           annotation (extent=[0, 40; 20, 20], rotation=180);
      
      Interfaces.JunctionVolume junctionVolume3(
        T_start=from_degC(50.0),
        V=V,
      redeclare package Medium = Modelica_Media.Water.IF97_ph) 
        annotation (extent=[-30, -50; -10, -70], rotation=180);
      Sources.FixedAmbient fixedAmbient4(
        p_ambient=from_bar(1.0),
        T_ambient=from_degC(10),
      redeclare package Medium = Modelica_Media.Water.IF97_ph) 
        annotation (extent=[-100, -70; -80, -50]);
      Components.ShortPipe shortPipe4(
        m_dot_nominal=10,
        frictionType=Modelica_Fluid.Types.FrictionTypes.ConstantLaminar,
      redeclare package Medium = Modelica_Media.Water.IF97_ph) 
        annotation (extent=[-60, -70; -40, -50]);
      Components.ShortPipe shortPipe5(
        m_dot_nominal=10,
        frictionType=Modelica_Fluid.Types.FrictionTypes.ConstantLaminar,
      redeclare package Medium = Modelica_Media.Water.IF97_ph) 
        annotation (extent=[0, -70; 20, -50]);
      Sources.FixedAmbient fixedAmbient5(
        T_ambient=from_degC(20),
        p_ambient=from_bar(0.95),
      redeclare package Medium = Modelica_Media.Water.IF97_ph) 
        annotation (extent=[60, -70; 40, -50]);
      Sources.FixedAmbient fixedAmbient6(
        p_ambient=from_bar(1.0),
        T_ambient=from_degC(20),
      redeclare package Medium = Modelica_Media.Water.IF97_ph) 
        annotation (extent=[-100, -30; -80, -10]);
      Components.ShortPipe shortPipe6(
        m_dot_nominal=10,
        frictionType=Modelica_Fluid.Types.FrictionTypes.ConstantLaminar,
      redeclare package Medium = Modelica_Media.Water.IF97_ph) 
        annotation (extent=[-60, -30; -40, -10]);
    equation 
      connect(fixedAmbient1.port, shortPipe1.port_a) 
        annotation (points=[-79, 30; -61, 30], style(color=69));
      connect(shortPipe3.port_b, fixedAmbient3.port) 
        annotation (points=[61, 30; 79, 30], style(color=69));
      connect(fixedAmbient2.port, shortPipe2.port_a) 
        annotation (points=[-79, 70; -61, 70], style(color=69));
      connect(shortPipe2.port_b, junctionVolume1.port) 
        annotation (points=[-39, 70; -20, 70; -20, 30], style(color=69));
      connect(shortPipe1.port_b, junctionVolume1.port) 
        annotation (points=[-39, 30; -20, 30], style(color=69));
      connect(junctionVolume2.port, shortPipe3.port_a) 
        annotation (points=[10, 30; 39, 30], style(color=69));
      connect(junctionVolume1.port, junctionVolume2.port) 
        annotation (points=[-20, 30; 10, 30], style(color=69));
      connect(fixedAmbient4.port, shortPipe4.port_a) 
        annotation (points=[-79, -60; -61, -60], style(color=69));
      connect(shortPipe5.port_b, fixedAmbient5.port) 
        annotation (points=[21, -60; 39, -60], style(color=69));
      connect(shortPipe5.port_a, junctionVolume3.port) 
        annotation (points=[-1, -60; -20, -60], style(color=69));
      connect(fixedAmbient6.port, shortPipe6.port_a) 
        annotation (points=[-79, -20; -61, -20], style(color=69));
      connect(shortPipe6.port_b, junctionVolume3.port) 
        annotation (points=[-39, -20; -20, -20; -20, -60], style(color=69));
      connect(shortPipe4.port_b, junctionVolume3.port) 
        annotation (points=[-39, -60; -20, -60], style(color=69));
    end TwoVolumesDetailedWater;
    
    model OneVolumeDetailedWater 
      "This example shows the difference of a JunctionVolume and of a FixedComponentVolume"
      
      
      import SI = Modelica.SIunits;
      parameter SI.Volume V=1.e-4 "Size of volume";
      
      import Modelica.SIunits.Conversions.*;
      
      extends Modelica.Icons.Example;
      annotation (Diagram, Documentation(info="<html>
<p>
This example shows two ways of merging flows. In the upper part
a volume in the connection port is modeled. At the initial time the
fluid in this volume has a higher temperature as the fluid
flowing form the ambient. Therefore, im shortPipe the temperature
will need some time to arrive at the mixing temperature.
</p>
<p>
In the lower part the volume in the connection port is neglected
and ideal mixing takes place. Since there is a steady state mass
flow and no (existing) delay mechanism is modeled, the mixing
is instantaneous. The temperature is the same as in the upper part.
</p>
</html>
"));
      Interfaces.JunctionVolume junctionVolume3(
        T_start=from_degC(50.0),
        V=V,
        redeclare package Medium = Modelica_Media.Water.IF97_ph) 
        annotation (extent=[-10, 0; 10, -20], rotation=180);
      Sources.FixedAmbient fixedAmbient4(
        p_ambient=from_bar(1.0),
        T_ambient=from_degC(10),
        redeclare package Medium = Modelica_Media.Water.IF97_ph) 
        annotation (extent=[-80, -20; -60, 0]);
      Components.ShortPipe shortPipe4(
        m_dot_nominal=10,
        frictionType=Modelica_Fluid.Types.FrictionTypes.ConstantLaminar,
        redeclare package Medium = Modelica_Media.Water.IF97_ph) 
        annotation (extent=[-40, -20; -20, 0]);
      Components.ShortPipe shortPipe5(
        m_dot_nominal=10,
        frictionType=Modelica_Fluid.Types.FrictionTypes.ConstantLaminar,
        redeclare package Medium = Modelica_Media.Water.IF97_ph) 
        annotation (extent=[20, -20; 40, 0]);
      Sources.FixedAmbient fixedAmbient5(
        T_ambient=from_degC(20),
        p_ambient=from_bar(0.95),
        redeclare package Medium = Modelica_Media.Water.IF97_ph) 
        annotation (extent=[80, -20; 60, 0]);
      Sources.FixedAmbient fixedAmbient6(
        p_ambient=from_bar(1.0),
        T_ambient=from_degC(20),
        redeclare package Medium = Modelica_Media.Water.IF97_ph) 
        annotation (extent=[-80, 20; -60, 40]);
      Components.ShortPipe shortPipe6(
        m_dot_nominal=10,
        frictionType=Modelica_Fluid.Types.FrictionTypes.ConstantLaminar,
        redeclare package Medium = Modelica_Media.Water.IF97_ph) 
        annotation (extent=[-40, 20; -20, 40]);
    equation 
      connect(fixedAmbient4.port, shortPipe4.port_a) 
        annotation (points=[-59, -10; -41, -10], style(color=69));
      connect(shortPipe5.port_b, fixedAmbient5.port) 
        annotation (points=[41, -10; 59, -10], style(color=69));
      connect(shortPipe5.port_a, junctionVolume3.port) 
        annotation (points=[19, -10; 0, -10], style(color=69));
      connect(fixedAmbient6.port, shortPipe6.port_a) 
        annotation (points=[-59, 30; -41, 30], style(color=69));
      connect(shortPipe6.port_b, junctionVolume3.port) 
        annotation (points=[-19, 30; 0, 30; 0, -10], style(color=69));
      connect(shortPipe4.port_b, junctionVolume3.port) 
        annotation (points=[-19, -10; 0, -10], style(color=69));
    end OneVolumeDetailedWater;
    
    model TwoVolumesApproximationWater 
      "This example shows the difference of a JunctionVolume and of a FixedComponentVolume"
      
      
      import SI = Modelica.SIunits;
      parameter SI.Volume V=1.e-4 "Size of volume";
      
      Interfaces.JunctionVolume junctionVolume1(
        T_start=from_degC(50.0),
        V=V/2,
        p_start=from_bar(100.0),
        initType=Modelica_Media.Interfaces.PartialMedium.Choices.Init.
            InitialStates, 
      redeclare package Medium = Modelica_Media.Water.IF97LocalApproximation_ph)
        annotation (extent=[-12,40; 8,20],     rotation=180);
      
      import Modelica.SIunits.Conversions.*;
      
      extends Modelica.Icons.Example;
      Sources.FixedAmbient fixedAmbient1(
        p_ambient=from_bar(100.0),
        T_ambient=from_degC(110), 
      redeclare package Medium = Modelica_Media.Water.IF97LocalApproximation_ph)
        annotation (extent=[-100, 20; -80, 40]);
      
      Components.ShortPipe shortPipe1(
        m_dot_nominal=10,
        frictionType=Modelica_Fluid.Types.FrictionTypes.ConstantLaminar,
        T_start=400.0, 
      redeclare package Medium = Modelica_Media.Water.IF97LocalApproximation_ph)
        annotation (extent=[-60, 20; -40, 40]);
      
      annotation (Diagram, Documentation(info="<html>
<p>
This example shows two ways of merging flows. In the upper part
a volume in the connection port is modeled. At the initial time the
fluid in this volume has a higher temperature as the fluid
flowing form the ambient. Therefore, im shortPipe the temperature
will need some time to arrive at the mixing temperature.
</p>
<p>
In the lower part the volume in the connection port is neglected
and ideal mixing takes place. Since there is a steady state mass
flow and no (existing) delay mechanism is modeled, the mixing
is instantaneous. The temperature is the same as in the upper part.
</p>
</html>
"));
      Components.ShortPipe shortPipe3(
        m_dot_nominal=10,
        frictionType=Modelica_Fluid.Types.FrictionTypes.ConstantLaminar,
        T_start=400, 
      redeclare package Medium = Modelica_Media.Water.IF97LocalApproximation_ph)
        annotation (extent=[40, 20; 60, 40]);
      
      Sources.FixedAmbient fixedAmbient3(
        p_ambient=from_bar(98),
        T_ambient=from_degC(120), 
      redeclare package Medium = Modelica_Media.Water.IF97LocalApproximation_ph)
        annotation (extent=[100, 20; 80, 40]);
      
    equation 
      connect(fixedAmbient1.port, shortPipe1.port_a) 
        annotation (points=[-79, 30; -61, 30], style(color=69));
      connect(shortPipe3.port_b, fixedAmbient3.port) 
        annotation (points=[61, 30; 79, 30], style(color=69));
      connect(shortPipe1.port_b, junctionVolume1.port) 
        annotation (points=[-39,30; -2,30],    style(color=69));
      connect(junctionVolume1.port, shortPipe3.port_a) annotation (points=[-2,
            30; 39,30], style(color=69, rgbcolor={0,127,255}));
    end TwoVolumesApproximationWater;
  end TestComponents;
  
  package Elementary 
    "Elementary examples to demonstrate various features of the MultiBody library"
    
    
    extends Modelica.Icons.Library;
    model SimpleMixing 
      "This example shows the difference of a JunctionVolume and of a FixedComponentVolume"
      
      
      Interfaces.JunctionVolume junctionVolume(
        V=1.e-4,
        T_start=from_degC(50.0),
        redeclare package Medium = Modelica_Media.Water.SimpleLiquidWater) 
        annotation (extent=[-10, 40; 10, 20], rotation=180);
      
      import Modelica.SIunits.Conversions.*;
      
      extends Modelica.Icons.Example;
      Sources.FixedAmbient fixedAmbient1(
        p_ambient=from_bar(1.0),
        T_ambient=from_degC(0),
        redeclare package Medium = Modelica_Media.Water.SimpleLiquidWater) 
        annotation (extent=[-80, 20; -60, 40]);
      Components.ShortPipe shortPipe1(
        m_dot_nominal=10,
        redeclare package Medium = Modelica_Media.Water.SimpleLiquidWater,
        frictionType=Modelica_Fluid.Types.FrictionTypes.ConstantLaminar) 
        annotation (extent=[-40, 20; -20, 40]);
      annotation (Diagram, Documentation(info="<html>
<p>
This example shows two ways of merging flows. In the upper part
a volume in the connection port is modeled. At the initial time the
fluid in this volume has a higher temperature as the fluid
flowing form the ambient. Therefore, im shortPipe the temperature
will need some time to arrive at the mixing temperature.
</p>
<p>
In the lower part the volume in the connection port is neglected
and ideal mixing takes place. Since there is a steady state mass
flow and no (existing) delay mechanism is modeled, the mixing
is instantaneous. The temperature is the same as in the upper part.
</p>
</html>
"));
      Components.ShortPipe shortPipe3(
        m_dot_nominal=10,
        redeclare package Medium = Modelica_Media.Water.SimpleLiquidWater,
        frictionType=Modelica_Fluid.Types.FrictionTypes.ConstantLaminar) 
        annotation (extent=[20, 20; 40, 40]);
      Sources.FixedAmbient fixedAmbient3(
        T_ambient=from_degC(20),
        p_ambient=from_bar(0.95),
        redeclare package Medium = Modelica_Media.Water.SimpleLiquidWater) 
        annotation (extent=[80, 20; 60, 40]);
      Sources.FixedAmbient fixedAmbient2(
        p_ambient=from_bar(1.0),
        redeclare package Medium = Modelica_Media.Water.SimpleLiquidWater,
        T_ambient=from_degC(10)) annotation (extent=[-80, 60; -60, 80]);
      Components.ShortPipe shortPipe2(
        m_dot_nominal=10,
        redeclare package Medium = Modelica_Media.Water.SimpleLiquidWater,
        frictionType=Modelica_Fluid.Types.FrictionTypes.ConstantLaminar) 
        annotation (extent=[-40, 60; -20, 80]);
      Sources.FixedAmbient fixedAmbient4(
        p_ambient=from_bar(1.0),
        redeclare package Medium = Modelica_Media.Water.SimpleLiquidWater,
        T_ambient=from_degC(1)) annotation (extent=[-80, -80; -60, -60]);
      Components.ShortPipe shortPipe4(
        m_dot_nominal=10,
        redeclare package Medium = Modelica_Media.Water.SimpleLiquidWater,
        frictionType=Modelica_Fluid.Types.FrictionTypes.ConstantLaminar) 
        annotation (extent=[-40, -80; -20, -60]);
      Components.ShortPipe shortPipe6(
        m_dot_nominal=10,
        redeclare package Medium = Modelica_Media.Water.SimpleLiquidWater,
        frictionType=Modelica_Fluid.Types.FrictionTypes.ConstantLaminar) 
        annotation (extent=[20, -80; 40, -60]);
      Sources.FixedAmbient fixedAmbient6(
        T_ambient=from_degC(20),
        p_ambient=from_bar(0.95),
        redeclare package Medium = Modelica_Media.Water.SimpleLiquidWater) 
        annotation (extent=[80, -80; 60, -60]);
      Sources.FixedAmbient fixedAmbient5(
        p_ambient=from_bar(1.0),
        redeclare package Medium = Modelica_Media.Water.SimpleLiquidWater,
        T_ambient=from_degC(10)) annotation (extent=[-80, -40; -60, -20]);
      Components.ShortPipe shortPipe5(
        m_dot_nominal=10,
        redeclare package Medium = Modelica_Media.Water.SimpleLiquidWater,
        frictionType=Modelica_Fluid.Types.FrictionTypes.ConstantLaminar) 
        annotation (extent=[-40, -40; -20, -20]);
    equation 
      connect(fixedAmbient1.port, shortPipe1.port_a) 
        annotation (points=[-59, 30; -41, 30], style(color=69));
      connect(shortPipe3.port_b, fixedAmbient3.port) 
        annotation (points=[41, 30; 59, 30], style(color=69));
      connect(shortPipe3.port_a, junctionVolume.port) 
        annotation (points=[19, 30; 0, 30], style(color=69));
      connect(fixedAmbient2.port, shortPipe2.port_a) 
        annotation (points=[-59, 70; -41, 70], style(color=69));
      connect(fixedAmbient4.port, shortPipe4.port_a) 
        annotation (points=[-59, -70; -41, -70], style(color=69));
      connect(shortPipe6.port_b, fixedAmbient6.port) 
        annotation (points=[41, -70; 59, -70], style(color=69));
      connect(fixedAmbient5.port, shortPipe5.port_a) 
        annotation (points=[-59, -30; -41, -30], style(color=69));
      connect(shortPipe5.port_b, shortPipe6.port_a) annotation (points=[-19, -30;
             0, -30; 0, -70; 19, -70], style(color=69));
      connect(shortPipe4.port_b, shortPipe6.port_a) 
        annotation (points=[-19, -70; 19, -70], style(color=69));
      connect(shortPipe2.port_b, junctionVolume.port) 
        annotation (points=[-19, 70; 0, 70; 0, 30], style(color=69));
      connect(shortPipe1.port_b, junctionVolume.port) 
        annotation (points=[-19, 30; 0, 30], style(color=69));
    end SimpleMixing;
    
  end Elementary;
  
  package Tanks "Examples with Tanks" 
    extends Icons.Library;
    
    model ThreeTanksOneLiquid 
      import Modelica.SIunits.Conversions.*;
      extends Modelica.Icons.Example;
      annotation (
        Diagram,
        Coordsys(grid=[1, 1], component=[20, 20]),
        experiment(StopTime=5),
        Documentation(info="<html>
<p>
This example demonstrates the mixing of a single substance flow
between three tanks with different temperature.
</p>
<p>
When comparing with Examples.Tanks.ThreeTanks, it is demonstrated
that it is easy in this case to switch between single and multiple
substance flow. The only differences between the examples ThreeTanks
and ThreeTanksOneLiquid are that the Medium is different and that
for the model ThreeTanksOneLiquid the default values are
used for the initial mass fractions.
</p>
</html>"));
      
      Components.Tank Tank1(
        area=1,
        T_start=from_degC(50),
        level_start=3,
        redeclare package Medium = Modelica_Media.Water.SimpleLiquidWater) 
        annotation (extent=[-90, 20; -70, 40]);
      Components.Tank Tank2(
        area=1,
        T_start=from_degC(100),
        level_start=1,
        redeclare package Medium = Modelica_Media.Water.SimpleLiquidWater) 
        annotation (extent=[-10, 20; 10, 40]);
      Components.ShortPipe shortPipe1(
        m_dot_nominal=2000,
        dp_nominal=from_bar(0.1),
        redeclare package Medium = Modelica_Media.Water.SimpleLiquidWater,
        frictionType=Modelica_Fluid.Types.FrictionTypes.ConstantLaminar) 
        annotation (extent=[-50, -30; -30, -10]);
      
      Components.Tank Tank3(
        area=1,
        T_start=from_degC(20),
        level_start=2,
        redeclare package Medium = Modelica_Media.Water.SimpleLiquidWater) 
        annotation (extent=[70, 20; 90, 40]);
      Components.ShortPipe shortPipe3(
        m_dot_nominal=2000,
        dp_nominal=from_bar(0.1),
        redeclare package Medium = Modelica_Media.Water.SimpleLiquidWater,
        frictionType=Modelica_Fluid.Types.FrictionTypes.ConstantLaminar) 
        annotation (extent=[30, -30; 50, -10]);
      Components.ShortPipe shortPipe2(
        m_dot_nominal=1000,
        dp_nominal=from_bar(0.1),
        redeclare package Medium = Modelica_Media.Water.SimpleLiquidWater,
        frictionType=Modelica_Fluid.Types.FrictionTypes.ConstantLaminar) 
        annotation (extent=[-10, -10; 10, 10], rotation=-90);
    equation 
      connect(Tank1.port, shortPipe1.port_a) 
        annotation (points=[-80, 19; -80, -20; -51, -20], style(color=69));
      connect(shortPipe3.port_b, Tank3.port) 
        annotation (points=[51, -20; 80, -20; 80, 19], style(color=69));
      connect(Tank2.port, shortPipe2.port_a) 
        annotation (points=[0, 19; 0, 11; -6.73533e-016, 11], style(color=69));
      connect(shortPipe1.port_b, shortPipe3.port_a) 
        annotation (points=[-29, -20; 29, -20], style(color=69));
      connect(shortPipe2.port_b, shortPipe3.port_a) annotation (points=[
            6.73533e-016, -11; 0, -11; 0, -20; 29, -20], style(color=69));
    end ThreeTanksOneLiquid;
    
    model ThreeTanksWithPortVolume 
      import Modelica.SIunits.Conversions.*;
      extends Modelica.Icons.Example;
      annotation (
        Diagram,
        Coordsys(grid=[1, 1], component=[20, 20]),
        experiment(StopTime=5),
        Documentation(info="<html>
<p>
This example demonstrates the mixing of a single substance flow
between three tanks with different temperatures. The difference to
example \"ThreeTanksOneLiquid\" is that the port where the three
pipes are connected together contains a volume now, in order that
the mixing of the pipe flows is modelled more realistically.
</p>
</html>"));
      
      Components.Tank Tank1(
        area=1,
        T_start=from_degC(50),
        level_start=3,
        redeclare package Medium = Modelica_Media.Water.SimpleLiquidWater) 
        annotation (extent=[-90, 20; -70, 40]);
      Components.Tank Tank2(
        area=1,
        T_start=from_degC(100),
        level_start=1,
        redeclare package Medium = Modelica_Media.Water.SimpleLiquidWater) 
        annotation (extent=[-10, 20; 10, 40]);
      Components.ShortPipe shortPipe1(
        m_dot_nominal=2000,
        dp_nominal=from_bar(0.1),
        redeclare package Medium = Modelica_Media.Water.SimpleLiquidWater,
        frictionType=Modelica_Fluid.Types.FrictionTypes.ConstantLaminar) 
        annotation (extent=[-50, -40; -30, -20]);
      
      Components.Tank Tank3(
        area=1,
        T_start=from_degC(20),
        level_start=2,
        redeclare package Medium = Modelica_Media.Water.SimpleLiquidWater) 
        annotation (extent=[70, 20; 90, 40]);
      Components.ShortPipe shortPipe3(
        m_dot_nominal=2000,
        dp_nominal=from_bar(0.1),
        redeclare package Medium = Modelica_Media.Water.SimpleLiquidWater,
        frictionType=Modelica_Fluid.Types.FrictionTypes.ConstantLaminar) 
        annotation (extent=[29, -40; 49, -20]);
      Components.ShortPipe shortPipe2(
        m_dot_nominal=1000,
        dp_nominal=from_bar(0.1),
        redeclare package Medium = Modelica_Media.Water.SimpleLiquidWater,
        frictionType=Modelica_Fluid.Types.FrictionTypes.ConstantLaminar) 
        annotation (extent=[-10, -10; 10, 10], rotation=-90);
      Interfaces.JunctionVolume junctionVolume(
        redeclare package Medium = Modelica_Media.Water.SimpleLiquidWater,
        V=1.e-4,
        T_start=from_degC(50.0),
        initType=Modelica_Media.Interfaces.PartialMedium.Choices.Init.
            InitialStates) annotation (extent=[-10, -40; 10, -20]);
    equation 
      connect(Tank1.port, shortPipe1.port_a) 
        annotation (points=[-80, 19; -80, -30; -51, -30], style(color=69));
      connect(shortPipe3.port_b, Tank3.port) 
        annotation (points=[50, -30; 80, -30; 80, 19], style(color=69));
      connect(Tank2.port, shortPipe2.port_a) 
        annotation (points=[0, 19; 0, 11; -6.73533e-016, 11], style(color=69));
      connect(shortPipe2.port_b, junctionVolume.port) annotation (points=[
            6.73533e-016, -11; 0, -11; 0, -30], style(color=69));
      connect(shortPipe1.port_b, junctionVolume.port) 
        annotation (points=[-29, -30; 0, -30], style(color=69));
      connect(shortPipe3.port_a, junctionVolume.port) 
        annotation (points=[28, -30; 0, -30], style(color=69));
    end ThreeTanksWithPortVolume;
    
    model ThreeTanksIF97 
      import Modelica.SIunits.Conversions.*;
      extends Modelica.Icons.Example;
      annotation (
        Diagram,
        Coordsys(grid=[1, 1], component=[20, 20]),
        experiment(StopTime=5),
        Documentation(info="<html>
<p>
This example demonstrates the mixing of a single substance flow
between three tanks with different temperature.
</p>
<p>
When comparing with Examples.Tanks.ThreeTanks, it is demonstrated
that it is easy in this case to switch between single and multiple
substance flow. The only differences between the examples ThreeTanks
and ThreeTanksOneLiquid are that the Medium is different and that
for the model ThreeTanksOneLiquid the default values are
used for the initial mass fractions.
</p>
</html>"));
      
      Components.Tank Tank1(
        area=1,
        T_start=from_degC(50),
        level_start=3,
        redeclare package Medium = Modelica_Media.Water.IF97_ph) 
        annotation (extent=[-90, 20; -70, 40]);
      Components.Tank Tank2(
        area=1,
        T_start=from_degC(100),
        level_start=1,
        redeclare package Medium = Modelica_Media.Water.IF97_ph) 
        annotation (extent=[-10, 20; 10, 40]);
      Components.ShortPipe shortPipe1(
        m_dot_nominal=2000,
        dp_nominal=from_bar(0.1),
        frictionType=Modelica_Fluid.Types.FrictionTypes.ConstantLaminar,
        redeclare package Medium = Modelica_Media.Water.IF97_ph) 
        annotation (extent=[-50, -30; -30, -10]);
      
      Components.Tank Tank3(
        area=1,
        T_start=from_degC(20),
        level_start=2,
        redeclare package Medium = Modelica_Media.Water.IF97_ph) 
        annotation (extent=[70, 20; 90, 40]);
      Components.ShortPipe shortPipe3(
        m_dot_nominal=2000,
        dp_nominal=from_bar(0.1),
        frictionType=Modelica_Fluid.Types.FrictionTypes.ConstantLaminar,
        redeclare package Medium = Modelica_Media.Water.IF97_ph) 
        annotation (extent=[30, -30; 50, -10]);
      Components.ShortPipe shortPipe2(
        m_dot_nominal=1000,
        dp_nominal=from_bar(0.1),
        frictionType=Modelica_Fluid.Types.FrictionTypes.ConstantLaminar,
        redeclare package Medium = Modelica_Media.Water.IF97_ph) 
        annotation (extent=[-10, -10; 10, 10], rotation=-90);
    equation 
      connect(Tank1.port, shortPipe1.port_a) 
        annotation (points=[-80, 19; -80, -20; -51, -20], style(color=69));
      connect(shortPipe3.port_b, Tank3.port) 
        annotation (points=[51, -20; 80, -20; 80, 19], style(color=69));
      connect(Tank2.port, shortPipe2.port_a) 
        annotation (points=[0, 19; 0, 11; -6.73533e-016, 11], style(color=69));
      connect(shortPipe1.port_b, shortPipe3.port_a) 
        annotation (points=[-29, -20; 29, -20], style(color=69));
      connect(shortPipe2.port_b, shortPipe3.port_a) annotation (points=[
            6.73533e-016, -11; 0, -11; 0, -20; 29, -20], style(color=69));
    end ThreeTanksIF97;
  end Tanks;
  
  package DrumBoiler 
    "Drum boiler example, see Franke, Rode, Krueger: On-line Optimization of Drum Boiler Startup, 3rd International Modelica Conference, Linkoping, 2003"
    
    
    model DrumBoilerSimulation "Simulate start-up of DrumBoiler" 
      extends Modelica.Icons.Example;
      Components.DrumBoiler drumBoiler annotation (extent=[-20, -40; 40, 20]);
      Modelica.Blocks.Sources.TimeTable q_F_Tab(table=[0, 0; 3600, 400; 7210,
            400]) annotation (extent=[-80, 2; -60, 22]);
      Modelica.Blocks.Sources.TimeTable Y_Valve_Tab(table=[0, 1; 3600, 1; 7210,
             1]) annotation (extent=[-80, -40; -60, -20]);
      annotation (
        Diagram,
        experiment(StopTime=7200),
        Documentation(info="<HTML>
<p>
Apply a ramp to fuel input and hold outlet valve open.
Simulate for 7200 seconds.
</p>
</HTML>"));
    equation 
      connect(q_F_Tab.outPort, drumBoiler.q_F) annotation (points=[-59, 12; -40,
             12; -40, -16; -25.7, -16], style(color=3));
      connect(Y_Valve_Tab.outPort, drumBoiler.Y_Valve) annotation (points=[-59,
             -30; -44, -30; -44, -34; -25.7, -34], style(
          color=3,
          fillColor=7,
          fillPattern=1));
    end DrumBoilerSimulation;
    
    package Components 
      
      model Evaporator 
        "Simple Evaporator with two states, see Astroem, Bell: Drum-boiler dynamics, Automatica 36, 2000, pp.363-378"
        
        
        import Modelica_Fluid.Interfaces.*;
        import Modelica.SIunits.Conversions.*;
        import SI = Modelica.SIunits;
        
        // property and interface declarations
        package Medium = WaterPhaseBoundaryIF97;
        Medium.BaseProperties medium_a(region=1, p=port_a.p) "Medium in port_a";
        Medium.BaseProperties medium_b(region=2, p=port_b.p) "Medium in port_b";
        FluidPort_a port_a(redeclare package Medium = Medium) 
          annotation (extent=[-120, -10; -100, 10]);
        FluidPort_b port_b(redeclare package Medium = Medium) 
          annotation (extent=[120, -10; 100, 10]);
        Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatPort 
          annotation (extent=[-10, -120; 10, -100]);
        Modelica.Blocks.Interfaces.OutPort V(redeclare type SignalType = 
              SI.Volume) "liquid volume (level)" 
          annotation (extent=[30, 100; 50, 120], rotation=90);
        Modelica.Blocks.Interfaces.OutPort sigma_D "Thermal stress in metal" 
          annotation (extent=[100, 40; 120, 60]);
        annotation (
          Coordsys(grid=[1, 1], component=[20, 20]),
          Diagram,
          Icon(
            Rectangle(extent=[-100, 59; 100, -61], style(
                color=0,
                gradient=2,
                fillColor=8)),
            Rectangle(extent=[-100, 34; 100, -36], style(
                color=69,
                gradient=2,
                fillColor=69)),
            Ellipse(extent=[18, 0; 48, -29], style(pattern=0, fillColor=7)),
            Ellipse(extent=[-1, 29; 29, 0], style(pattern=0, fillColor=7)),
            Ellipse(extent=[48, 34; 78, 5], style(pattern=0, fillColor=7)),
            Ellipse(extent=[-31, 1; -1, -28], style(pattern=0, fillColor=7)),
            Ellipse(extent=[47, 14; 77, -15], style(pattern=0, fillColor=7)),
            Ellipse(extent=[-72, 25; -42, -4], style(pattern=0, fillColor=7)),
            Ellipse(extent=[71, 0; 101, -29], style(pattern=0, fillColor=7)),
            Ellipse(extent=[74, 14; 104, -15], style(pattern=0, fillColor=7)),
            Ellipse(extent=[71, 29; 101, 0], style(pattern=0, fillColor=7)),
            Text(
              extent=[-120, 117; 116, 51],
              string="%name",
              style(gradient=2, fillColor=69)),
            Line(points=[0, -60; 0, -100], style(color=42)),
            Line(points=[40, 99; 40, 60])));
        
        // public parameters
        parameter SI.Mass m_D=300e3 "mass of surrounding drum metal";
        parameter SI.SpecificHeatCapacity cp_D=500 
          "specific heat capacity of drum metal";
        parameter SI.Volume V_t=100 "total volume inside drum";
        parameter SI.Pressure p_start=from_bar(1) "initial pressure";
        parameter SI.Volume V_start=67 "initial liquid volume";
        
      protected 
        SI.Pressure p(start=p_start, stateSelect=StateSelect.prefer) 
          "pressure inside drum boiler";
        SI.Volume V_v "volume of vapour phase";
        SI.Volume V_l(start=V_start, stateSelect=StateSelect.prefer) 
          "volumes of liquid phase";
        SI.SpecificEnthalpy h_v=medium_b.h "specific enthalpy of vapour";
        SI.SpecificEnthalpy h_l=medium_a.h "specific enthalpy of liquid";
        SI.Density rho_v=medium_b.d "density in vapour phase";
        SI.Density rho_l=medium_a.d "density in liquid phase";
        SI.Mass m "total mass of drum boiler";
        SI.Energy U "internal energy";
        SI.Temperature T_D=heatPort.T "temperature of drum";
        SI.HeatFlowRate q_F=heatPort.Q_dot "heat flow rate from furnace";
        SI.SpecificEnthalpy h_W=port_a.h "feed water enthalpy";
        SI.SpecificEnthalpy h_S=medium_b.h "steam enthalpy";
        SI.MassFlowRate qm_W=port_a.m_dot "feed water mass flow rate";
        SI.MassFlowRate qm_S=port_b.m_dot "steam mass flow rate";
      equation 
        
        // balance equations  
        m = rho_v*V_v + rho_l*V_l + m_D;
        U = rho_v*V_v*h_v + rho_l*V_l*h_l - p*V_t + m_D*cp_D*T_D;
        der(m) = qm_W + qm_S;
        der(U) = q_F + qm_W*h_W + qm_S*h_S;
        T_D = medium_a.T;
        // ideal heat transfer between metal and water
        V_t = V_l + V_v;
        
        // pressure and specific total enthalpies at ports
        port_a.p = p;
        port_b.p = p;
        port_b.H_dot = semiLinear(port_b.m_dot, port_b.h, h_v);
        port_a.H_dot = semiLinear(port_a.m_dot, port_a.h, h_l);
        
        // thermal stress
        sigma_D.signal[1] = 60*der(T_D);
        
        // liquid level 
        V.signal[1] = V_l;
      end Evaporator;
      
      model Valve 
        "Simple controlled valve with linear pressure drop coefficient" 
        
        import SI = Modelica.SIunits;
        import Modelica.SIunits.Conversions.*;
        import Modelica_Fluid.*;
        extends Interfaces.PartialTwoPortTransport;
        SI.Pressure dp "Pressure loss due to friction";
        Real residue=port_a.p - port_b.p - dp 
          "momentum balance (may be modified)";
        
        parameter Real k=1e-5 "linear valve coefficient";
        
        Modelica.Blocks.Interfaces.InPort Y "Valve position" 
          annotation (extent=[-10, -80; 10, -60], rotation=90);
        
        annotation (Icon(
            Text(
              extent=[-126, -76; 130, -110],
              style(color=0),
              string=""),
            Text(
              extent=[-120, 130; 116, 64],
              string="%name",
              style(gradient=2, fillColor=69)),
            Line(points=[-60, -50; -60, 50; 60, -50; 60, 50; -60, -50], style(
                  color=0, thickness=2)),
            Line(points=[-60, 0; -100, 0], style(color=69)),
            Line(points=[60, 0; 100, 0], style(color=69)),
            Line(points=[0, 0; 0, -72])));
      equation 
        residue = 0;
        port_a.m_dot = Y.signal[1]*k*dp;
      end Valve;
      
      model MassFlowSource_h 
        "Ideal pump that produces a mass flow rate from a large reservoir defined by input signal"
        
        
        import Modelica_Fluid.*;
        import Modelica_Fluid.Interfaces.*;
        import Modelica_Media.Interfaces.*;
        
        import SI = Modelica.SIunits;
        import Modelica.SIunits.Conversions.*;
        
        FluidPort_b port(redeclare model Medium = Medium) 
          annotation (extent=[100, -10; 120, 10], rotation=0);
        replaceable package Medium = PartialMedium extends PartialMedium 
          "Medium in the component"             annotation (choicesAllMatching=true);
        
        parameter SI.SpecificEnthalpy h_ambient=5e5 "Ambient enthalphy";
        parameter SI.Pressure p_start=from_bar(1.0) 
          "|Initialization|| Initial pressure";
        annotation (
          Coordsys(
            extent=[-100, -100; 100, 100],
            grid=[2, 2],
            component=[20, 20]),
          Icon(
            Rectangle(extent=[20, 60; 100, -60], style(
                color=0,
                gradient=2,
                fillColor=8)),
            Rectangle(extent=[38, 40; 100, -40], style(
                color=69,
                gradient=2,
                fillColor=69)),
            Ellipse(extent=[-100, 80; 60, -80], style(fillColor=7)),
            Polygon(points=[-60, 70; 60, 0; -60, -68; -60, 70], style(color=73,
                   fillColor=73)),
            Text(
              extent=[-54, 32; 16, -30],
              string="m'",
              style(color=41, fillColor=41)),
            Text(extent=[-142, 142; 156, 88], string="%name"),
            Text(
              extent=[-126, -86; 146, -116],
              style(color=0),
              string="%h_ambient")),
          Window(
            x=0.45,
            y=0.01,
            width=0.44,
            height=0.65),
          Diagram);
        Modelica.Blocks.Interfaces.InPort m_dot(redeclare type SignalType = 
              SI.MassFlowRate) 
          "Mass flow rate from an infinite reservoir in to the port as signal" 
          annotation (extent=[-140, -20; -100, 20]);
      equation 
        port.m_dot = -noEvent(max(m_dot.signal[1], 0));
        port.H_dot = semiLinear(port.m_dot, port.h, h_ambient);
      end MassFlowSource_h;
      
      model DrumBoiler 
        "Complete drum boiler model, including evaporator and supplementary components"
        
        
        import Modelica.SIunits.Conversions.*;
        
        Evaporator evaporator annotation (extent=[-50, -20; -30, 0]);
        annotation (
          uses(Modelica_Fluid(version="0.72")),
          Diagram,
          Coordsys(
            extent=[-100, -100; 100, 100],
            grid=[1, 1],
            component=[20, 20]),
          Icon(
            Rectangle(extent=[-100, 100; 100, -100], style(fillColor=7)),
            Text(
              extent=[-151, 165; 138, 102],
              style(fillColor=7, fillPattern=1),
              string="%name"),
            Text(
              extent=[-79, 67; 67, 21],
              string="drum",
              style(
                color=0,
                fillColor=7,
                fillPattern=1)),
            Text(
              extent=[-90, -14; 88, -64],
              string="boiler",
              style(
                color=0,
                fillColor=7,
                fillPattern=1))));
        Modelica.Thermal.HeatTransfer.PrescribedHeatFlow furnace 
          annotation (extent=[-50, -50; -30, -30], rotation=90);
        Modelica.Blocks.Interfaces.InPort q_F 
          annotation (extent=[-139, 0; -99, -40]);
        Modelica.Blocks.Interfaces.InPort Y_Valve 
          annotation (extent=[-139, -100; -99, -60]);
        Valve valve(redeclare package Medium = Modelica_Media.Water.IF97_ph, k=
              1.5e-5) annotation (extent=[44, -20; 64, 0]);
        Modelica_Fluid.Sources.FixedAmbient sink(redeclare package Medium = 
              Modelica_Media.Water.IF97_ph) 
          annotation (extent=[80, -20; 100, 0], rotation=180);
        Modelica_Fluid.Sensors.MassFlowSensor massFlowSensor(redeclare package 
            Medium = Modelica_Media.Water.IF97_ph) 
          annotation (extent=[10, -20; 30, 0]);
        Modelica_Fluid.Sensors.TemperatureSensor temperatureSensor(redeclare 
            package Medium = Modelica_Media.Water.IF97_ph) 
          annotation (extent=[10, 60; 30, 80]);
        Modelica_Fluid.Sensors.PressureSensor pressureSensor(redeclare package 
            Medium = Modelica_Media.Water.IF97_ph) 
          annotation (extent=[10, 20; 30, 40]);
        Modelica.Blocks.Continuous.PI controller(k={10}, T={120}) 
          annotation (extent=[-60, 30; -80, 50]);
        MassFlowSource_h pump(redeclare package Medium = 
              Modelica_Media.Water.IF97_ph) 
          annotation (extent=[-80, -20; -60, 0]);
        Modelica.Blocks.Math.Feedback feedback 
          annotation (extent=[-26, 30; -46, 50]);
        Modelica.Blocks.Sources.Constant levelSetPoint(k={67}) 
          annotation (extent=[-46, 60; -26, 80]);
      protected 
        Modelica.Blocks.Interfaces.OutPort T_S 
          annotation (extent=[34, 66; 42, 74]);
        Modelica.Blocks.Interfaces.OutPort p_S 
          annotation (extent=[34, 26; 42, 34]);
        Modelica.Blocks.Interfaces.OutPort qm_S 
          annotation (extent=[34, -34; 42, -26], rotation=0);
        Modelica.Blocks.Interfaces.OutPort sigma_D 
          annotation (extent=[-24, 6; -16, 14]);
        Modelica.Blocks.Interfaces.OutPort V_l 
          annotation (extent=[-24, 24; -16, 32]);
      public 
        Modelica.Blocks.Math.Gain MW2W(k={1e6}) 
          annotation (extent=[-70, -69; -50, -50]);
      equation 
        connect(furnace.port, evaporator.heatPort) 
          annotation (points=[-40, -30; -40, -21], style(color=42));
        connect(Y_Valve, valve.Y) 
          annotation (points=[-119, -80; 54, -80; 54, -17], style(color=3));
        connect(evaporator.port_b, temperatureSensor.port) annotation (points=[
              -29, -10; -2, -10; -2, 50; 20, 50; 20, 59], style(color=69));
        connect(evaporator.port_b, pressureSensor.port) annotation (points=[-29,
               -10; -2, -10; -2, 10; 20, 10; 20, 19], style(color=69));
        connect(evaporator.port_b, massFlowSensor.port_a) 
          annotation (points=[-29, -10; 9, -10], style(color=69));
        connect(massFlowSensor.port_b, valve.port_a) 
          annotation (points=[31, -10; 43, -10], style(color=69));
        connect(valve.port_b, sink.port) annotation (points=[65, -10; 72, -10;
              72, -10; 79, -10], style(color=69));
        connect(pump.port, evaporator.port_a) 
          annotation (points=[-59, -10; -51, -10], style(color=69));
        connect(controller.inPort, feedback.outPort) 
          annotation (points=[-58, 40; -45, 40], style(color=3));
        connect(feedback.inPort2, evaporator.V) 
          annotation (points=[-36, 32; -36, 1], style(color=3));
        connect(levelSetPoint.outPort, feedback.inPort1) annotation (points=[-25,
               70; -20, 70; -20, 40; -28, 40], style(color=3));
        connect(pressureSensor.signal, p_S) 
          annotation (points=[30, 30; 38, 30], style(color=3));
        connect(temperatureSensor.signal, T_S) 
          annotation (points=[30, 70; 38, 70], style(color=3));
        connect(massFlowSensor.signal, qm_S) 
          annotation (points=[20, -20; 20, -30; 38, -30], style(color=3));
        connect(evaporator.sigma_D, sigma_D) annotation (points=[-29, -5; -26,
              -5; -26, 10; -20, 10], style(color=3));
        connect(evaporator.V, V_l) 
          annotation (points=[-36, 1; -36, 28; -20, 28], style(color=3));
        connect(controller.outPort, pump.m_dot) annotation (points=[-81, 40; -90,
               40; -90, -10; -82, -10], style(color=3));
        connect(q_F, MW2W.inPort) annotation (points=[-119, -20; -90, -20; -90,
               -59.5; -72, -59.5], style(color=3));
        connect(MW2W.outPort, furnace.Q_dot) annotation (points=[-49, -59.5; -40,
               -59.5; -40, -50], style(color=3));
      end DrumBoiler;
      
      package WaterPhaseBoundaryIF97 
        "Physical properties for water at phase boundary at boiling and dew curves"
        
        
        extends Modelica_Media.Interfaces.PartialMedium(
          mediumName="WaterIF97",
          substanceNames=fill("", 0),
          incompressible=false,
          reducedX=true,
          MassFlowRate(quantity="MassFlowRate.WaterIF97"));
        
        redeclare model extends BaseProperties 
          
          parameter Integer region=0 "specify region 1 (liquid) or 2 (vapour)";
        equation 
          
          assert(region == 1 or region == 2,
            "WaterPhaseBoundaryIF97 medium model only valid for regions 1 and 2");
          T = Modelica_Media.Water.IF97.BaseIF97.Basic.tsat(p);
          if region == 1 then
            d = Modelica_Media.Water.IF97.BaseIF97.Regions.rhol_p(p);
            h = Modelica_Media.Water.IF97.BaseIF97.Regions.hl_p(p);
          else
            d = Modelica_Media.Water.IF97.BaseIF97.Regions.rhov_p(p);
            h = Modelica_Media.Water.IF97.BaseIF97.Regions.hv_p(p);
          end if;
          u = h - p/d;
        end BaseProperties;
      end WaterPhaseBoundaryIF97;
    end Components;
  end DrumBoiler;
  
end Examples;


package Interfaces 
  "Interfaces for steady state and unsteady, mixed-phase, multi-substance, incompressible and compressible flow without momentum."
  
  
  annotation (Documentation(info="<HTML>
<p>
This library provides connector definitions for the following
type of flow systems:
</p>
<ul>
<li> One-dimensional flow of a <b>single substance</b>
     or of a <b>mixture of substances</b> with optional <b>multiple phases</b>.</li>
<li> <b>Incompressible</b> and <b>compressible</b> medium.</li>
<li> <b>Steady state</b> and <b>unsteady</b> flow.</li>
<li> <b>Momentum</b> of the flow is <b>neglected</b>.</li>
<li> The <b>kinetic energy</b> in the flow is <b>neglected</b><br>
     (this is usually a good approximation, if |v| &lt; 0.1 Ma, i.e.,
     if the flow speed is less than about 10 % of the speed of sound in
     the medium).</li>
</ul>
<p>
Additionally, the medium model interfaces are defined with package
<b>PartialMedium</b>. The definitions are made in such a way that
it is not possible to connect connectors of different media together.
</p>
<dl>
<dt><b>Main Author:</b>
<dd>Hilding Elmqvist, Dynasim</dl>
<p><b>Release Notes:</b></p>
<ul>
<li><i>Nov. 6, 2002</i>
       by Hilding Elmqvist: first version.</li>
<li><i>Nov. 11, 2002</i>
       by Hilding Elmqvist, Martin Otter: improved version.</li>
<li><i>Nov. 20-21, 2002</i>
       by Hilding Elmqvist, Mike Tiller, Allan Watson, John Batteh, Chuck Newman,
       Jonas Eborn: Improved at the 32nd Modelica Design Meeting.
<li><i>Aug. 11, 2002</i>
       by Martin Otter: Improved according to discussion with Hilding
       Elmqvist and Hubertus Tummescheit.<br>
       The PortVicinity model is manually
       expanded in the base models.<br>
       The Volume used for components is renamed
       PartialComponentVolume.<br>
       A new volume model \"Fluid.Components.PortVolume\"
       introduced that has the medium properties of the port to which it is
       connected.<br>
       Fluid.Interfaces.PartialTwoPortTransport is a component
       for elementary two port transport elements, whereas PartialTwoPort
       is a component for a container component.</li>
</li>
</ul>
</HTML>"));
  
  extends Modelica.Icons.Library;
  import SI = Modelica.SIunits;
  
  connector FluidPort 
    "Interface for quasi one-dimensional fluid flow in a piping network (incompressible or compressible, one or more phases, one or more substances)"
    
    
    replaceable package Medium = Modelica_Media.Interfaces.PartialMedium 
      "Medium model" annotation (choicesAllMatching=true);
    
    Medium.AbsolutePressure p "Pressure in the connection point";
    flow Medium.MassFlowRate m_dot 
      "Mass flow rate from the connection point into the component";
    
    Medium.SpecificEnthalpy h 
      "Specific mixture enthalpy in the connection point";
    flow Medium.EnthalpyFlowRate H_dot 
      "Enthalpy flow rate into the component (if m_dot > 0, H_dot = m_dot*h)";
    
    Medium.MassFraction X[Medium.nX](quantity=Medium.substanceNames) 
      "Independent mixture mass fractions m_i/m in the connection point";
    flow Medium.MassFlowRate mX_dot[Medium.nX](quantity=Medium.substanceNames) 
      "Mass flow rates of the independent substances from the connection point into the component (if m_dot > 0, mX_dot = m_dot*X)";
    
  end FluidPort;
  
  connector FluidPort_a "Fluid connector with filled icon" 
    extends FluidPort;
    annotation (Diagram(Rectangle(extent=[-100, 100; 100, -100], style(color=69,
               fillColor=69)), Text(extent=[-88, 206; 112, 112], string="%name")),
         Icon(Rectangle(extent=[-100, 100; 100, -100], style(color=69,
              fillColor=69)), Text(
          extent=[-126, 160; 130, 104],
          string="%name",
          style(
            color=0,
            fillColor=69,
            fillPattern=1))));
  end FluidPort_a;
  
  connector FluidPort_b "Fluid connector with outlined icon" 
    extends FluidPort;
    annotation (Diagram(Rectangle(extent=[-100, 100; 100, -100], style(color=69,
               fillColor=7)), Text(extent=[-88, 192; 112, 98], string="%name")),
         Icon(Rectangle(extent=[-100, 100; 100, -100], style(color=69,
              fillColor=7)), Text(
          extent=[-126, 160; 130, 104],
          string="%name",
          style(
            color=0,
            fillColor=69,
            fillPattern=1))));
  end FluidPort_b;
  
  partial model PartialInit 
    "Define Medium model and parameter menu to initialize medium in component that has states"
    
    
    import Modelica.SIunits.Conversions.*;
    replaceable package Medium = PackageMedium extends 
      Modelica_Media.Interfaces.PartialMedium "Medium in the component" 
                  annotation (choicesAllMatching=true);
    
    parameter Medium.Choices.Init.Temp initType=Medium.Choices.Init.NoInit 
      "|Initialization||Type of initialization" annotation (Evaluate=true);
    parameter Boolean init_p=true 
      "|Initialization|Initialization of mass balance| = true, if p_start is used, otherwise d_start (true is required for incompressible medium)";
    parameter Medium.AbsolutePressure p_start=101325 
      "|Initialization|Initialization of mass balance| Start value of pressure p, if init_p = true";
    parameter Medium.Density d_start=1 
      "|Initialization|Initialization of mass balance| Start value of density d, if init_p = false";
    parameter Boolean init_T=true 
      "|Initialization|Initialization of energy balance| = true, if T_start is used, otherwise h_start";
    parameter Medium.Temperature T_start=from_degC(20) 
      "|Initialization|Initialization of energy balance| Start value of temperature T, if init_T = true";
    parameter Medium.SpecificEnthalpy h_start=1.e4 
      "|Initialization|Initialization of energy balance| Start value of specific enthalpy h, if init_T = false";
    parameter Medium.MassFraction X_start[:]=zeros(Medium.nX) 
      "|Initialization|Initialization of mass fractions (only for multi-substance fluids)| Start values of independent mass fractions X";
  end PartialInit;
  
  partial model PartialInitAlgebraic 
    "Define Medium model and parameter menu to initialize medium for algebraic equations (e.g. ShortPipe)"
    
    
    import Modelica.SIunits.Conversions.*;
    
    replaceable package Medium = PackageMedium extends 
      Modelica_Media.Interfaces.PartialMedium "Medium in the component" 
                  annotation (choicesAllMatching=true);
    
    parameter Boolean init_p=true 
      "|Initialization|Initialization of mass balance| = true, if p_start is used, otherwise d_start (true is required for incompressible medium)";
    parameter Medium.AbsolutePressure p_start=101325 
      "|Initialization|Initialization of mass balance| Start value of pressure p, if init_p = true";
    parameter Medium.Density d_start=1.e3 
      "|Initialization|Initialization of mass balance| Start value of density d, if init_p = false";
    parameter Boolean init_T=true 
      "|Initialization|Initialization of energy balance| = true, if T_start is used, otherwise h_start";
    parameter Medium.Temperature T_start=from_degC(20) 
      "|Initialization|Initialization of energy balance| Start value of temperature T, if init_T = true";
    parameter Medium.SpecificEnthalpy h_start=1.e4 
      "|Initialization|Initialization of energy balance| Start value of specific enthalpy h, if init_T = false";
    parameter Medium.MassFraction X_start[:]=zeros(Medium.nX) 
      "|Initialization|Initialization of mass fractions (only for multi-substance fluids)| Start values of independent mass fractions X";
  end PartialInitAlgebraic;
  
  partial model PartialSource 
    "Partial component source with one fluid connector" 
    
    replaceable package Medium = PackageMedium extends 
      Modelica_Media.Interfaces.PartialMedium 
      "Medium model within the source (= different from the port medium, if port.m_dot > 0)"
                   annotation (choicesAllMatching=true);
    
    FluidPort_b port(redeclare package Medium = Medium) 
      annotation (extent=[100, -10; 120, 10], rotation=0);
    
    Medium.BaseProperties medium "Medium in the source";
  equation 
    port.p = medium.p;
    
      /* Handle reverse and zero flow (for details, see Fluid.Interfaces.semiLinear) */
    port.H_dot = semiLinear(port.m_dot, port.h, medium.h);
    port.mX_dot = semiLinear(port.m_dot, port.X, medium.X);
    annotation (Documentation(info="<html>
<p>
Partial component to model the <b>volume interface</b> of a <b>source</b>
component, such as a mass flow source. The essential
features are:
</p>
<ul>
<li> The pressure in the connection port (= port.p) is identical to the
     pressure in the volume (= medium.p).</li>
<li> The enthalpy flow rate (= port.H_dot) and the mass flow rates of the
     substances (= port.mX_dot) depend on the direction of the mass flow rate,
     according to the semiLinear(..) equations.</li>
</ul>
</html>"));
  end PartialSource;
  
  partial model PartialTwoPortTransport 
    "Partial element transporting fluid between two ports without storing mass or energy"
    
    
    import SI = Modelica.SIunits;
    import Cv = Modelica.SIunits.Conversions;
    
    extends PartialInitAlgebraic;
    
    FluidPort_a port_a(redeclare package Medium = Medium) 
      annotation (extent=[-120, -10; -100, 10]);
    FluidPort_b port_b(redeclare package Medium = Medium) 
      annotation (extent=[120, -10; 100, 10]);
    Medium.BaseProperties medium_a(
      final init_p=init_p,
      final p_start=p_start,
      final d_start=d_start,
      final init_T=init_T,
      final T_start=T_start,
      final h_start=h_start,
      final X_start=X_start) "Medium properties in port_a";
    Medium.BaseProperties medium_b(
      final init_p=init_p,
      final p_start=p_start,
      final d_start=d_start,
      final init_T=init_T,
      final T_start=T_start,
      final h_start=h_start,
      final X_start=X_start) "Medium properties in port_b";
    Medium.MassFlowRate m_dot 
      "Mass flow rate from port_a to port_b (m_dot > 0 is design flow direction)";
    
    annotation (
      Coordsys(grid=[1, 1], component=[20, 20]),
      Diagram,
      Documentation(info="<html>
<p>
This component transports fluid between its two ports, without
storing mass or energy. Reversal and zero mass flow rate is taken
care off, for details see definition of built-in operator semiLinear().
When using this partial component,
the momentum equation has to be added by specifying a relationship
between the pressure drop \"dp = port_a.p - port_b.p\" and the
mass flow rate \"port_a.m_dot\".
</p>
</html>"));
  equation 
    /* Compute medium variables from port variables */
    medium_a.p = port_a.p;
    medium_a.h = port_a.h;
    medium_a.X = port_a.X;
    medium_b.p = port_b.p;
    medium_b.h = port_b.h;
    medium_b.X = port_b.X;
    
    /* Handle reverse and zero flow */
    port_a.H_dot = semiLinear(port_a.m_dot, port_a.h, port_b.h);
    port_a.mX_dot = semiLinear(port_a.m_dot, port_a.X, port_b.X);
    
    /* Energy, mass and substance mass balance */
    port_a.H_dot + port_b.H_dot = 0;
    port_a.m_dot + port_b.m_dot = 0;
    port_a.mX_dot + port_b.mX_dot = zeros(Medium.nX);
    
    // Design direction of mass flow rate
    m_dot = port_a.m_dot;
  end PartialTwoPortTransport;
  
  partial model PartialOnePort 
    "Partial fluid component with one fluid connector (to be used as container of other components)"
    
    
    extends PartialInit;
    
    FluidPort_a port(redeclare package Medium = Medium) 
      annotation (extent=[-10, -120; 10, -100], rotation=90);
    annotation (
      Icon,
      Coordsys(grid=[1, 1], component=[20, 20]),
      Diagram,
      Documentation(info="<html>
<p>
This partial component should be used for new <b>container</b> models
with <b>two fluid ports</b> that contain
connections of elementary fluid component models (such as Fluid.Components.Pipe or
Fluid.Components.PortVolume).
</p>
<p>
Note, whenever a new elementary component
is implemented that is defined directly by equations, reversal and zero mass
flow rate has to be taken care off. This is most easily performed by using
the built-in operator semiLinear(). When connecting elementary
components together, the flow reversal is already handeled in these components.
</p>
</html>"));
  end PartialOnePort;
  
  partial model PartialTwoPort 
    "Partial fluid component with two fluid connectors (to be used as container of other components)"
    
    
    extends PartialInit;
    
    FluidPort_a port_a(redeclare package Medium = Medium) 
      annotation (extent=[-120, -10; -100, 10]);
    FluidPort_b port_b(redeclare package Medium = Medium) 
      annotation (extent=[120, -10; 100, 10]);
    annotation (
      Coordsys(grid=[1, 1], component=[20, 20]),
      Diagram,
      Documentation(info="<html>
<p>
This partial component should be used for new <b>container</b> models with <b>two fluid ports</b>
that contain connections of elementary fluid component models (such as Fluid.Components.Pipe or
Fluid.Components.PortVolume).
</p>
<p>
Note, whenever a new elementary component
is implemented that is defined directly by equations, reversal and zero mass
flow rate has to be taken care off. This is most easily performed by using
the submodel Modelica_Fluid.Interfaces.semiLinear. When connecting elementary
components together, the flow reversal is already handeled in these components.
</p>
</html>"));
  end PartialTwoPort;
  
  model ThermalAdaptor 
    "Adaptor between a HeatPort (heat transfer port) and a FluidPort" 
    
    replaceable package Medium = PackageMedium "Medium model" annotation (
        choicesAllMatching=true);
    Medium.BaseProperties medium "Medium used in the this model";
    
    FluidPort_a fluidPort(redeclare package Medium = Medium) 
      annotation (extent=[100, -10; 120, 10], rotation=0);
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatPort 
      annotation (extent=[-120, -10; -100, 10]);
    annotation (
      Coordsys(
        extent=[-100, -100; 100, 100],
        grid=[2, 2],
        component=[20, 20]),
      Icon(
        Text(extent=[-134, 134; 134, 72], string="%name"),
        Rectangle(extent=[0, 60; 60, -60], style(
            color=0,
            gradient=1,
            fillColor=9)),
        Rectangle(extent=[-20, 60; 0, -60], style(color=0, fillColor=42)),
        Rectangle(extent=[-100, 6; -20, -6], style(color=42, fillColor=42)),
        Rectangle(extent=[60, 60; 100, -60], style(
            color=69,
            fillColor=69,
            fillPattern=1))),
      Diagram,
      Documentation(info="<html>
<p>
This component is an adaptor between a connector from the Modelica.Thermal.HeatTransfer
library and a connector from the Fluid library. This adaptor is used
to model the heat transfer of a part into a fluid by connecting the
heat transfer components via this adaptor to a connection point of
the fluid.
</p>
</html>"));
  equation 
    
      // Intensive quantities of the fluidPort are used to compute medium properties
    medium.p = fluidPort.p;
    medium.h = fluidPort.h;
    medium.X = fluidPort.X;
    
    // No mass flow from the heatPort to the fluidPort
    fluidPort.m_dot = 0;
    fluidPort.mX_dot = zeros(Medium.nX);
    
    // Energy balance between the two ports
    heatPort.Q_dot + fluidPort.H_dot = 0;
    
    // Boundary condition
    heatPort.T = medium.T;
  end ThermalAdaptor;
  
  model JunctionVolume 
    "Fixed volume associated with a port by the finite volume method (the medium properties of the volume are the ones of the port)"
    
    
    import SI = Modelica.SIunits;
    import Cv = Modelica.SIunits.Conversions;
    
    extends PartialInit;
    
    FluidPort_b port(redeclare package Medium = Medium) 
      annotation (extent=[-10, -10; 10, 10], rotation=0);
    
    parameter SI.Volume V=1e-6 "Fixed size of junction volume";
    
    Medium.BaseProperties medium(
      preferedMediumStates=true,
      final initType=initType,
      final init_p=init_p,
      final p_start=p_start,
      final d_start=d_start,
      final init_T=init_T,
      final T_start=T_start,
      final h_start=h_start,
      final X_start=X_start);
    
    SI.Energy U "Internal energy of port volume";
    Real m(quantity=Medium.mediumName, unit="kg") "Mass of junction volume";
    Real mX[Medium.nX](quantity=Medium.substanceNames, each unit="kg") 
      "Independent substance masses of junction volume";
  equation 
    medium.p = port.p;
    medium.h = port.h;
    medium.X = port.X;
    m = V*medium.d;
    U = m*medium.u;
    mX = m*medium.X;
    der(m) = port.m_dot;
    der(U) = port.H_dot;
    der(mX) = port.mX_dot;
    annotation (Icon(
        Ellipse(extent=[-100, 100; 100, -100], style(
            color=0,
            gradient=3,
            fillColor=71)),
        Text(extent=[-144, 178; 146, 116], string="%name"),
        Text(
          extent=[-130, -108; 144, -150],
          style(color=0),
          string="V=%V")), Documentation(info="<html>
<p>
This component models the <b>volume</b> of <b>fixed size</b> that is
associated with the <b>fluid port</b> to which it is connected.
This means that all medium properties inside the volume, are identical
to the port medium properties. In particular, the specific enthalpy
inside the volume (=h) is always identical to the specific enthalpy
in the port (port.h = h). Usually, this model is used when
discretizing a component according to the finite volume method into
volumes in internal ports that only store energy and mass and into
transport elements that just transport energy, mass and momentum
between the internal ports without storing these quantities during the
transport.
</p>
</html>"));
    
  end JunctionVolume;
  
  function checkAmbient "Check whether ambient definition is correct" 
    
    extends Modelica.Icons.Function;
    input String mediumName;
    input Boolean incompressible;
    input Boolean define_p;
    input Boolean reducedX;
    input Integer nX;
    input Real X_ambient[:];
  algorithm 
    assert(not incompressible or incompressible and define_p, "
Wrong value of parameter define_p (= false) in ambient source component:
The selected medium \"" + mediumName + "\" is incompressible.
Therefore, an ambient density cannot be defined and
define_p = true is required.
");
    
    for i in 1:nX loop
      assert(X_ambient[i] >= 0.0, "
Wrong ambient mass fractions in medium \"" + mediumName + "\":
The ambient value X_ambient(" + integerString(i) + ") = " + realString(
        X_ambient[i]) + "
is negative. It must be positive.
");
    end for;
    
    assert(reducedX or not reducedX and nX > 0 and abs(sum(X_ambient) - 1.0) <
      1.e-10, "
Wrong ambient mass fractions in medium \"" + mediumName + "\":
This medium requires that the ambient mass fractions X_ambient
sum up to 1. However, sum(X_ambient) = " + realString(sum(X_ambient)) + ".
");
    
    assert(not reducedX or reducedX and sum(X_ambient) < 1 + 1.e-10, "
Wrong ambient mass fractions in medium \"" + mediumName + "\":
This medium requires that the sum of the ambient mass fractions X_ambient
is at most 1. However, sum(X_ambient) = " + realString(sum(X_ambient)) + ".
");
  end checkAmbient;
end Interfaces;
=======
>>>>>>> 1.2
end Modelica_Fluid;
