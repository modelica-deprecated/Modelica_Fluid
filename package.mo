package Modelica_Fluid "1-dimensional thermo-fluid flow in pipe networks using the Modelica.Media media description (requires package Modelica 2.2.2 development)"
annotation (
  version="1.0 Beta 3",
  versionDate="2007-06-05",
  preferedView="info",
  Settings(NewStateSelection=true),
  uses(Modelica(version="2.2.1")),
  classOrder={"UsersGuide","Examples","Ambient", "ControlValves","Flowmachines","HeatExchangers","Junctions",
      "Volumes", "Pipes", "PressureLosses", "Pumps", "Sensors", "Sources", "Thermal", "*"},
  Documentation(info="<html>
<p>
Library <b>Modelica_Fluid</b> is a <b>free</b> Modelica package providing
components describing
<b>1-dimensional thermo-fluid flow</b> in networks of pipes. A unique feature is that the
component equations and the media models are decoupled.
All components are implemented such that they can be used for
media from the Modelica.Media library. This means especially that an
incompressible or compressible medium, a single or a multiple
substance medium with one or more phases might be used.
The goal is to include 
the Modelica_Fluid library in the Modelica standard library as Modelica.Fluid.
</p>
 
<p>
This is version <b>1.0 Beta 3</b> of the Modelica_Fluid library.
We expect that the structure and the components of the library do 
not change or only marginally change for the 1.0 release version.
For the 1.0 version the documentation will be improved,
a few components might be slightly changed and the simulation
will be made more stable. Please, read the section
<a href=\"Modelica:Modelica_Fluid.UsersGuide.KnownLimitations\">Known limitations</a>
in the Users Guide before using this library.
</p>
 
<p>
A typical example model of the Modelica_Fluid library
is shown in the next figure (drum boiler):
</p>
<p align=\"center\">
<img src=\"../Images/UsersGuide/DrumBoiler.png\">
</p>
<p>
An example of a tank system that is controlled by a control system
and where some of the components have built-in diagram animation
is shown in the next figure:
</p>
<p align=\"center\">
<img src=\"../Images/Examples/ControlledTanks1.png\">
</p>
<p>
The following parts are useful, when newly starting with this library:
</p>
<ul>
<li> <a href=\"Modelica:Modelica_Fluid.UsersGuide\">Modelica_Fluid.UsersGuide</a>.</li>
<li> <a href=\"Modelica:Modelica_Fluid.UsersGuide.ReleaseNotes\">Modelica_Fluid.UsersGuide.ReleaseNotes</a>
     summarizes the changes of the library releases.</li>
<li> <a href=\"Modelica:Modelica_Fluid.Examples\">Modelica_Fluid.Examples</a>
     contains examples that demonstrate the usage of this library.</li>
</ul>
 
<p>
Note, Modelica_Fluid does <b>not</b> work with
version 2.2 of the Modelica standard library. 
The reason is that some additional functions have been
added to Modelica.Media in 2.2.1 that are accessed in Modelica_Fluid.
</p>
 
<p><b>Copyright &copy; 2002-2006, Modelica Association.</b></p>
<p><i>
This Modelica package is <b>free</b> software; it can be redistributed and/or modified
under the terms of the <b>Modelica license</b>, see the license conditions
and the accompanying <b>disclaimer</b> in the documentation of package
Modelica in file \"Modelica/package.mo\".
</i></p>
 
 
</html>"),
    conversion(from(version="0.795", script=
            "../ConvertFromModelica_Fluid_0.795.mos")));


  extends Modelica.Icons.Library;
  import SI = Modelica.SIunits;


package UsersGuide "Users Guide" 
  
  annotation (DocumentationClass=true, Documentation(info="<HTML>
<h3><font color=\"#008000\" size=5>Users guide of package Modelica_Fluid</font></h3>
<p>
Library <b>Modelica_Fluid</b> is a <b>free</b> Modelica package providing
components describing
<b>1-dimensional thermo-fluid flow</b> in networks of pipes. A unique feature is that the
component equations and the media models are decoupled.
All components are implemented such that they can be used for
media from the Modelica.Media library. This means especially that an
incompressible or compressible medium, a single or a multiple
substance medium with one or more phases might be used.
The goal is to include 
the Modelica_Fluid library in the Modelica standard library as Modelica.Fluid.
</p>
</HTML>"));
  
class KnownLimitations "Known limitations" 
    
    annotation (Documentation(info="<html>
<h3><font color=\"#008000\" size=5>Known limitations</font></h3>

<p>
This first public release of the Modelica_Fluid library is
still a Beta-Version. The goal is to improve it and we are
interested in your feedback, i.e., bug reports and
improvement suggestions. Please, send your emails to
the fluid development group at <b>Modelica-design@Modelica.org</b>.
</p>

<p>
The Modelica_Fluid library has quite ambitious goals, especially,
that every component can be connected in an arbitrary way,
e.g., pipes can be flipped, without influencing the generated
code, and that medium descriptions with different independent
variables can be used (e.g., \"p,T\" or \"p,h\" as medium variables).
This is a new approach in the fluid modeling area and as with
every new approach, it takes some time that everything works
as expected. We are aware of the following limitations of the current
version:
</p>

<ul>
<li> The medium has to be defined in <b>every</b> component
     separately. The user would like to define the medium
     at one location and then, the definition should propagate
     through the connections automatically. With our current
     general approach, this is <b>not</b> possible to formulate it
     in Modelica. There are discussions in the Modelica development
     group to extend Modelica with new featurres so that this
     becomes possible.<br>
     In Dymola, there is some type of tool support, to make this
     a bit better: Select all fluid components, right click with the
     mouse and select \"parameters\". Then select the desired medium
     and this medium will be used in all selected components.
    </li>

<li> When clicking on the <b>Medium</b> parameter, a very long list
     of media is displayed. We would like to have better control for
     the end user what is shown or at least display it hierarchically.
     It is not yet clear how this should be improved, but we work on it.
     </li>

<li> It might be that when connecting components together you get quickly
     large non-linear systems of equations and initialization might fail
     due to inappopriate start values. This will be improved.
     You can reduce the number of non-linear systems of equations (and thereby
     make the simulation more stable), by placing an instance of model
     \"Modelica_Fluid.Pipes.BaseClasses.PortVolume\" in every connection
     node (this means that a volume is present in every connection node).
     </li>
</ul>
</html>"));
end KnownLimitations;
  
  class Overview "Overview" 
    
    annotation (Documentation(info="<HTML>
<h3><font color=\"#008000\" size=5>Overview</font></h3>
<p>
The Modelica_Fluid library provides basic interfaces and 
components to model 1-dim. thermo-fluid flow in networks of pipes.
It is not the intention that this library covers all
application cases because the fluid flow area is too large and
because for special applications it is possible to implement
libraries with simpler component interfaces.
Instead, the goal is that the Modelica_Fluid library provides
a <b>reasonable set of components</b> and that it <b>demonstrates</b>
how to implement components of a fluid flow library in Modelica,
in particular to cope with difficult issues such as connector
design, reversing flow and initialization. It is planned to
include more components in the future. User proposals are
welcome.
</p>
<p>
This library has the following main features:
</p>
<ul>
<li> The connectors Modelica_Fluid.Interfaces.FluidPort_a/_b are designed
     for one-dimensional flow of a <b>single substance</b>
     or of a <b>mixture of substances</b> with optional <b>multiple phases</b>.
     All media models from Modelica.Media can be utilized when
     connecting components. For one substance media, the additional arrays for 
     multiple
     substance media have zero dimension and are therefore removed
     from the code during translation. The general connector definition
     therefore does not introduce an overhead for special cases.<br>&nbsp;</li>
<li> All the components of the Modelica_Fluid library are designed
     that they can be utilized for all media models from
     Modelica.Media if this is posssible. For example, all media can
     be utilized for the Modelica_Fluid.Sensors/Sources components.
     For some components only special media are possible, since additional
     functionality is required. For example, 
     Modelica_Fluid.Components.Evaporator requires a two phase medium
     (extending from Modelica.Media.Interfaces.PartialTwoPhaseMedium).
     <br>&nbsp;</li>
<li> In order to simplify the initialization in the components,
     there is the restriction that only media models are supported
     that have T, (p,T), (p,h), (T,X), (p,T,X) or (p,h,X) as
     independent variables. Other media models would be possible,
     e.g., with (T,d) as independent variables. However, this requires
     to rewrite the code for the component initialization.
     (Note, T is temperature, p is pressure, d is density,
     h is specific enthalpy, and X is a mass fraction vector).
     <br>&nbsp;</li>
<li> All components work for <b>incompressible</b> and <b>compressible</b> media. 
     This is implemented by a small change in the initialization of a 
     component, if the medium is incrompressible. Otherwise, the equations
     of the components are not influenced by this property.<br>&nbsp;</li>
<li> All components allow fluid flow in both directions, i.e., 
     <b>reversing flow</b> is supported.<br>&nbsp;</li> 
<li> 2 or more components can be connected together. The effect can
     be interpreted as introducing an infinitesimal small control
     volume in the connecting point with <b>ideal mixing</b>. In particular,
     the mass- and energy balance is automatically fulfilled.
     The momentum balance is only fulfilled if two components are
     connected in a straight line with the same pipe diameter.
     If more detailed
     models are needed for a connection, e.g., if mixing losses
     shall be taken into account, an appropriate model has to be 
     explicitly used in the connection point.<br>&nbsp;</li>
<li> There is no restriction how components can be connected
     together (besides from special cases).
     E.g., PortVolumes can be directly connnected
     to each other, as well as pressure drop components.
     Such types of connections might introduce additional
     linear or non-linear systems of equations.<br>&nbsp;</li>
     
</ul>
</HTML>
"));
  equation 
    
  end Overview;
  
  class GettingStarted "Getting started" 
    
    annotation (Documentation(info="<html>
<h3><font color=\"#008000\" size=5>Getting started</font></h3>
<p>
An example will be included here.
</p>
</html>
"));
  equation 
    
  end GettingStarted;
  
  class ComponentDefinition "Component definition" 
    
    annotation (Documentation(info="<html>
<h3><font color=\"#008000\" size=5>Component definition</font></h3>
<p>
In this section it is described how the components
of the Modelica_Fluid library are implemented.
If you would like to introduce new components either in
Modelica_Fluid or your own library, you should be aware
of the issues discussed in this section.
</p>
<p>
This section is partly based on the following paper:
</p>
<dl>
<dt> Elmqvist H., Tummescheit H., and Otter M.:</dt>
<dd> <b>Object-Oriented Modeling of Thermo-Fluid Systems</b>.
     Modelica 2003 Conference, Link&ouml;ping, Sweden, 
     pp. 269-286, Nov. 3-4, 2003.
     Download from:
     <a href=\"http://www.modelica.org/Conference2003/papers/h40_Elmqvist_fluid.pdf\">http://www.modelica.org/Conference2003/papers/h40_Elmqvist_fluid.pdf</a>
     </dd>
</dl>
</html>
"));
    
  class FluidConnectors "Fluid connectors" 
      
    annotation (Documentation(info="<html>
<h3><font color=\"#008000\" size=5>Fluid connectors</font></h3>
<p>
In this section the design of the fluid connectors is
explained. A major design goal was that components can be arbitrarily
connected and that the important balance equations are automatically
fulfilled when 2 or more components are connected together at
one point as shown in the next figure:
</p>
<p align=\"center\">
<img src=\"../Images/UsersGuide/MixingConnections.png\">
</p>
<p>
In such a case the balance equations define <b>ideal mixing</b>,
i.e., the connection point has the mixing temperature if the
fluids from the three components would be ideally mixed in 
an infinitely small time period. If more realistic modelling
is desired that takes into account mixing losses, an explicit
model has to be used in the connection point.
</p>
<h4><font color=\"#008000\">Single substance media</font></h4>
<p>
For a single substance medium, the connector definition in
Modelica_Fluid.Interfaces.FluidPort reduces to
</p>
<pre>
  <b>connector</b> FluidPort 
     <b>replaceable package</b> Medium = Modelica.Media.Interfaces.PartialMedium;
  
     Medium.AbsolutePressure  p;
     <b>flow</b> Medium.MassFlowRate m_flow;
  
     Medium.SpecificEnthalpy      h;
     <b>flow</b> Medium.EnthalpyFlowRate H_flow \"if m_flow &gt; 0, H_flow = m_flow*h\"
  <b>end</b> FluidPort;
</pre>
<p>
The first statement defines the Medium flowing through the connector.
In a medium, medium specific types such as \"Medium.AbsolutePressure\" 
are defined that contain medium specific values for the min, max and
nominal attributes. Furthermore, Medium.MassFlowRate is defined as:
</p>
<pre>
   <b>type</b> MassFlowRate = Modelica.SIunits.MassFlowRate(
                                    quantity=\"MassFlowRate.\" + mediumName, ...);
</pre>
<p>
A Modelica translator will check that the quantity and unit attributes 
of connected interfaces are identical. Therefore, an error occurs,
if connected FluidPorts do not have a medium with the same medium name.
</p>
<p>
As in all Modelica libraries, some requirements must be
fulfilled for a component, in order that the connection equations
generated by a Modelica translator from a \"connect(..)\" statement
lead to the balance equations in the particular area.
For fluid libraries, we have to first analyse the balance
equations present in a component. For one-dimensional flow
along the coordinate \"x\", the following partial differential
equations hold 
</p>
<table border=1 cellspacing=0 cellpadding=2>
  <tr><td> Mass balance</td>
      <td> <img src=\"../Images/UsersGuide/massBalance.png\"></td>
  </tr>
  <tr><td> Momentum balance</td>
      <td> <img src=\"../Images/UsersGuide/momentumBalance.png\"></td>
  </tr>
  <tr><td> Energy balance 1</td>
      <td> <img src=\"../Images/UsersGuide/energyBalance1.png\"></td>
  </tr>
  <tr><td> Pipe friction</td>
      <td> <img src=\"../Images/UsersGuide/pipeFriction.png\"></td>
  </tr>
  <tr><td></td>
      <td>x: independent spatial coordinate (flow is along coordinate x)<br>
          t: time<br>
          v(x,t): mean velocity<br>
          p(x,t): mean pressure<br>
          T(x,t): mean temperature<br>
          &rho;(x,t): mean density<br> 
          u(x,t): specific internal energy<br>
          z(x): height over ground<br>
          A(x): area perpendicular to direction x<br>
          g: gravity constant
          f: Fanning friction factor
          S: circumference
  </tr>
</table>
<p>
An alternative energy balance can be derived by multiplying
the momentum balance with \"v\" and substracting it
from the energy balance 1 above. This results in
the \"energy balance 2\":
</p>
<table border=1 cellspacing=0 cellpadding=2>
  <tr><td> Energy balance 2</td>
      <td> <img src=\"../Images/UsersGuide/energyBalance2.png\"></td>
  </tr>
</table>
<p>
This equation is much simpler because all terms depending
on the velocity v are removed, especially the kinetic
energy term. This means that the pure \"mechanical\" part of
the energy balance was removed since fulfilled identically
via the (mechanical) momentum equation and that the \"thermal\" part
of the energy balance remains. <b>All components</b> of the
Fluid library use the <b>energy balance 2</b> equation.
</p>
<p>
Assume that 3 components are connected together at one infinitesimal
small control volume and that thermal conduction
is neglected. Under the assumption of ideal mixing,
the intensive quantities at the ports of the components and in the
small control volume are identical. Furthermore, 
for the small control volume the 
mass and energy balance equations from above reduce
to the following simple equations (note, that neither
mass nor energy is stored in the volume and that
m_flow1 is the mass flow rate and
H_flow1 is the enthalpy flow rate into component 1):
</p>
<table border=1 cellspacing=0 cellpadding=2>
  <tr><td> Intensive quantities</td>
      <td> p1=p2=p3; h1=h2=h3; T1=T2=T3; etc.</td>
  <tr><td> Mass balance</td>
      <td> 0 = m_flow1 + m_flow2 + m_flow3</td>
  </tr>
  <tr><td> Energy balance 2</td>
      <td> 0 = H_flow1 + H_flow2 + H_flow3</td>
  </tr>
</table>
<p>
As can be seen these are exactly also the equations
generated by a Modelica translator from the connector
definition above: Non-flow variables are identical, i.e.,
p1=p2=p3, h1=h2=h3 and the flow variables sum up to zero.
Since the other intensive quantities, such as T or u, are
a function of two independent variables such as p and h, 
it follows that T1=T2=T3, u1=u2=u3 etc.
</p>
<p> 
A connector should have only the minimal number of variables to 
describe the interface, otherwise there will be connection
restrictions in certain cases. Therefore, in the connector
no redundant variables are present, e.g., the temperature T 
is not present because it can be computed from the connector
variables pressure p and specific enthalpy h (as will be
described in subsequent sections, this will not lead
to an increased computational effort, if appropriate
support by the Modelica translator is available).
</p>
<p>
The momentum equation in three dimensions reduces
to the following vector equation for the small
control volume:
</p>
<table border=1 cellspacing=0 cellpadding=2>
  <tr><td> 3D momentum balance</td>
      <td> 0 = m_flow1*<b>v1</b> + 
               m_flow2*<b>v2</b> + 
               m_flow3*<b>v2</b></td>
  </tr>
</table>
<p>
where <b>v1</b>, <b>v2</b>, <b>v3</b> are the
velocity vectors at the ports of the 3 components.
This equation is only fulfilled in certain cases.
For example, if two pipes are connected
together along a straight line with same pipe
areas, then <b>v1</b> = <b>v2</b> and the momentum
balance can be written as:
</p>
<pre>
    0 = <b>v1</b>*(m_flow1 + m_flow2)
      = 0
</pre>
<p>
because the term in paranthesis is the mass balance.
</p>
<p>
In the general case, the momentum equation is not fulfilled.
In several applications, it is a useful simplification to
neglect the momentum balance for a connecting point. In such
a case, 3 and more components can be directly connected.
In other cases, e.g., gas dynamics, the momentum balance
is essential and cannot be neglected. Then, a model
has to be used in the connection point and 3 and more
components cannot be directly connected together.
</p>
<p>
To summarize, all Fluid components shall be implemented
with the <b>energy balance 2</b> form, and the mass and energy
balance are fulfilled for the infinitesimal small control
volume in a connection point. The momentum balance is only
fulfilled in certain cases.  
If this is not justified, an own component has to be defined that 
models the connection, including the desired form of the momentum balance.
</p>
<h4><font color=\"#008000\">Multiple substance media</font></h4>
<p>
xxx
</p>
</html>
"));
  end FluidConnectors;
    
  class UpstreamDiscretization "Upstream discretization" 
      
    annotation (Documentation(info="<html>
<h3><font color=\"#008000\" size=5>Upstream discretization</font></h3>
<p>
When implementing a Fluid component, the difficult arises that
the value of intensive quantities (such as p, T, &rho;) 
shall be accessed from the
<b>upstream</b> volume. For example, if the fluid flows
from volume A to volume B, then the intensive quantities of
volume B have no influence on the fluid between the
two volumes. On the other hand, if the flow direction is reversed,
the intensive quantities of volume A have no influence
on the fluid between the two volumes.
</p>
<p>
In the Modelica_Fluid library, such a situation is handeled
with the following code fragment 
(from Interfaces.PartialTwoPortTransport,
simplified for a single substance fluid):
</p>
<pre>    <b>replaceable package</b> Medium =  
                   Modelica.Media.Interfaces.PartialMedium 
                   <b>annotation</b>(choicesAllMatching = <b>true</b>);
  
    Interfaces.FluidPort_a port_a(<b>redeclare package</b> Medium = Medium); 
    Interfaces.FluidPort_b port_b(<b>redeclare package</b> Medium = Medium); 
    Medium.BaseProperties medium_a \"Medium properties in port_a\";
    Medium.BaseProperties medium_b \"Medium properties in port_b\";
  <b>equation</b> 
    // Properties in the ports
    port_a.p = medium_a.p;
    port_a.h = medium_a.h;
    port_b.p = medium_b.p;
    port_b.h = medium_b.h;
  
    // Handle reverse and zero flow
    port_a.H_flow = <b>semiLinear</b>(port_a.m_flow, port_a.h, port_b.h);
  
    // Energy and mass balance; here: 
    port_a.H_flow + port_b.H_flow = 0;
    port_a.m_flow + port_b.m_flow = 0;
      ...
</pre>
<p>
The medium models medium_a and medium_b are associated with
the Fluid connectors port_a and port_b. The enthalpy
flow rate in port_a, port_a.H_flow, is in principal computed
with an if clause:
</p>
<pre>  port_a.H_flow = <b>if</b> port_a.m_flow &gt; 0 <b>then</b> port_a.m_flow*port_a.h
                                       <b>else</b> port_a.m_flow*port_b.h
                = <b>semiLinear</b>(port_a.m_flow, port_a.h,  port_b.h);
</pre>
<p>
However, instead of using this if-clause, the corresponding
built-in Modelica operator <b>semiLinear</b>() is used.
The main reason is that this operator will allow a Modelica
translator certain symbolic transformations that lead
to a more robust numerical computation. Note, the necessary
symbolic manipulation is only possible, if exactly the above
structure is used, i.e., one variable on the left of the
equal sign and the semiLinear() operator directly on the right
side of the equal sign. For example, the following construction
will not work and should therefore not be used:
</p>
<pre>   // wrong usage of semiLinear() operator:
   port_a.H_flow = cp*<b>semiLinear</b>(port_a.m_flow, medium_a.T, medium_b.T);
</pre>
<p>
If the above component is connected between two volumes, i.e.,
the independent medium variables in port_a and port_b are
usually states, then port_a.h and port_b.h are either states
(i.e., known quantities in the model) or are computed from
states. In either case they are \"known\". In such a situation, 
the above if-clause represented by the \"semiLinear\" operator
is uncritical, because it depends only on known variables and
can be directly computed. At zero mass flow rate,
a change from one branch of the if-clause to the other branch 
takes place. Since port_a.m_flow is zero, port_a.H_flow
changes <b>continuously</b> (only the derivative of the equation
is discontinuous at this point).
</p>
<p>
According to the Modelica language specification 2.2, an 
event shall be triggered at zero mass flow rate, because the
Boolean condition \"port_a.m_flow > 0\" changes its sign.
It turned out that this approach leads to problems in some
situations (unnecessary chattering around zero mass flow rate). 
Therefore, it was decided that for the next release of
the language, it is undefined whether an event is generated or
not. The recommendation is to not generate an event.
In Dymola, an event is not generated for Dymola version &ge; 5.3d+. 
The exact definition of the semiLinear() operator in the next
version of the Modelica language specification is given 
<a href=\"Modelica:Modelica_Fluid.UsersGuide.ComponentDefinition.SemiLinearDefinition\">here</a>.
</p>
<p>
If two components C1 and C2 of the above type are connected together
(e.g., two pressure drop pipe components), the following equations
are present (the semiLinear operators are replaced by the corresponding
if clauses and for simplicity, abbreviations of the form
H_flow1 = C1.port_a.H_flow are used; h is the specific enthalpy
in the connection point, h1 is the specific enthalpy in component C1
and h2 is the specific enthalpy in component C2):
</p>
<pre>  H_flow1 = <b>if</b> m_flow1 &ge; 0 <b>then</b> m_flow1*h <b>else</b> m_flow1*h1;
  H_flow2 = <b>if</b> m_flow2 &ge; 0 <b>then</b> m_flow2*h <b>else</b> m_flow1*h2;
        0 = m_flow1 + m_flow2;
        0 = H_flow1 + H_flow2;
</pre>
<p>
This equation system can be transformed to
</p>
<pre>   m_flow1 &ge; 0, m_flow2 &lt; 0:  0 = m_flow1*(h  - h2);  -&gt; h = h2
   m_flow1 &lt; 0, m_flow2 &ge; 0:  0 = m_flow1*(h1 - h );  -&gt; h = h1
</pre>
<p>
For m_flow1 = m_flow2 = 0, it follows from the original equations that
H_flow1 = H_flow2 = 0. This means that \"h\" is no longer appearing in an
equation and can be selected arbitrarily, i.e., there is an infinite
number of solutions. With the transformation rules for the 
semiLinear() operator, a Modelica translator can transform the 
original equations (if expressed with the semiLinear() operator)
into:
</p>
<pre>   h = <b>if</b> m_flow1 &ge; 0 <b>then</b> h2 <b>else</b> h1;
</pre>
<p>
In this case, the ambiguity is explicitly solved by (arbitrarily)
defining h=h2, in case m_flow1 = 0. Note, this is performed
during translation and therefore the ambiguity is no longer
present in the generated code. In reality, no ambiguity is present,
but the intrinsic thermal quantities at zero mass flow rate
are determined by thermal conduction. Since thermal conduction
along the fluid flow is neglected in Modelica_Fluid, this ambiguity
is present.
</p>
<p>
If 3 or more components are connected together, it is no longer
possible to resolve the ambiguity during translation. 
Instead, the systems of equations can be transformed
into the following scalar equation, where h is the unknown
specific enthalpy in the connection point:
</p>
<pre>  f1(m_flow1,m_flow2,m_flow3)*h = f2(m_flow1,m_flow2,m_flow3,h1,h2,h3); 
</pre>
<p>
If m_flow1 = m_flow2 = m_flow3 = 0, this equation degenerates into
</p>
<pre>  0*h = 0
</pre>
<p>
and therefore there is again an infinite number of solutions for h.
A Modelica simulation environment might resolve this ambiguity
by using the solution of h from the last accepted step.
</p>
  
</html>
"));
  end UpstreamDiscretization;
    
  class PropertyPropagation "Property propagation" 
      
    annotation (Documentation(info="<html>
<h3><font color=\"#008000\" size=5>Property propagation</font></h3>
<p>
As explained in section
<a href=\"Modelica:Modelica_Fluid.UsersGuide.ComponentDefinition.FluidConnectors\">Fluid connectors</a>,
it is possible to 
connect components together in a nearly arbitrary fashion,
because every connection fulfills automatically the
balance equations. This approach has, however, one drawback:
If two components are connected together, then the medium
variables on both sides of the connector are identical.
However, due to the connector, only the two equations
</p>
<pre>
   p1 = p2;
   h1 = h2;
</pre>
<p>
are present. Assume, that p,T are the independent medium variables
and that the medium properties are computed at one side of the
connections. This means, the following equations are basically
present:
</p>
<pre>
    h1 = h(p1,T1);
    h2 = h(p2,T2);
    p1 = p2;
    h1 = h2;
</pre>
<p>
These equations can be solved in the following way:
</p>
<pre>
    h1 := h(p1,T1)
    p2 := p1;
    h2 := h1;
    0  := h2 - h(p2,T2);   // non-linear system of equations for T2
</pre>
<p>
This means that T2 is computed by solving a non-linear system
of equations. If h1 and h2 are provided as Modelica functions,
a Modelica translator, such as Dymola, can replace
this non-linear system of equations by the equation:
</p>
<pre>
   T2 := T1;
</pre>
<p>
because after alias substition there are two function calls
</p>
<pre>
    h1 := h(p1,T1);
    h1 := h(p1,T2);
</pre>
<p>
Since the left hand side of the function call and the first
argument are the same, the second arguments T1 and T2 must also be 
identical and therefore T2 := T1. This type of analysis seems
to be only possible, if the specific enthalpy is defined as a function
of the independent medium variables. Due to this property, all
media in the Modelica.Media library define
the specific enthalpy always as a function and therefore by
appropriate tool support (as, e.g., in Dymola) no non-linear
system of equations appears and in the generated code
propagation of medium properties over a connector does not
lead to an unnecessary overhead.
</p>
</html>
"));
  end PropertyPropagation;
    
  class RegularizingCharacteristics "Regularizing characteristics" 
      
    annotation (Documentation(info="<html>
<h3><font color=\"#008000\" size=5>Regularizing characteristics</font></h3>
<p>
Pressure drop equations and other fluid characteristics are usually
computed by <b>semi-empirical</b> equations. Unfortunately, the developers
of semi-empirical equations nearly never take into account that the
equation might be used in a simulation program. As a consequence, these
semi-empirical equations can nearly never be used blindly but must
be slightly modified or adapted in order that obvious
simulation problems are avoided. Below, examples are given to
demonstrate what problems occur and how to regularize the characteristics:  
</p>
<h4><font color=\"#008000\">Square root function</font></h4>
<p>
In several empirical formulae, expressions of the following form
are present, e.g., for turbulent flow in a pipe:
</p>
<pre>   y = <b>if</b> x &lt; 0 <b>then</b> -<b>sqrt</b>( <b>abs</b>(x) ) <b>else</b> <b>sqrt</b>(x)
</pre>
<p>
A plot of this characteristic is shown in the next figure:
</p>
<p align=\"center\">
<img src=\"../Images/UsersGuide/sqrt.png\">
</p>
<p>
The difficulty with this function is that the derivative at x=0 is infinity.
In reality, such a function does not exist. E.g., for pipe flow,
the flow becomes laminar for small velocities and therefore around zero the
sqrt() function is replaced by a linear function. Since the laminar region is
usually of not much practical interest, the above approximation is used.
</p>
<p>
The direct implementation above does not work in Modelica, because
an event is generated when x &lt; 0 changes sign. In order to detect
this event, an event iteration takes place. During the event iteration,
the active if-branche is not changed. For example, assume that x is positive
(= \"else\" branch) and shall become negative. During the event iteration
x is slightly negative and the else branch, i.e., sqrt(x), is evaluated.
Since this results in an imaginary number, an error occurs.
It would be possible to fix this, by using the <b>noEvent</b>() operator
to explicitly switch of an event:
</p>
<pre>   y = <b>if</b> <b>noEvent</b>(x &lt; 0) <b>then</b> -<b>sqrt</b>( <b>abs</b>(x) ) <b>else</b> <b>sqrt</b>(x)
</pre>
<p>
Still, it is highly likely that good integrators will not work well
around x=0, because they will recognize that the derivative changes very
sharply and will reduce the step size drastically.
</p>
<p>
There are several solutions around this problem: Around x=0, the sqrt() function
can be replaced by a polynomial of 3rd order which is determined in such a way
that it smoothly touches the sqrt() function, i.e., the whole function is continuous
and continuously differentiable. In the Modelica_Fluid library, implementations of
such critical functions are provided in sublibrary Modelica_Fluid.Utilities.
The above sqrt() type function is computed by function <b>Utilities.regRoot</b>().
This function is defined as:
</p>
<pre>     y := x/(x*x+delta*delta)^0.25;
</pre>
<p>
where \"delta\" is the size of the small region around zero where the
sqrt() function is approximated by another function. The plot of the
function above is practically identical to the one of the original function.
However, it has a finite derivative at x=0 and is differentiable upto
any order. With the default value of delta=0.01, the difference between 
the function above and regRoot(x) is 16% around x=0.01, 0.25% around x=0.1 
and 0.0025% around x=1.
</p>
</html>
"));
  end RegularizingCharacteristics;
    
  class WallFriction "Wall friction" 
      
    annotation (Documentation(info="<html>
<h3><font color=\"#008000\" size=5>Wall friction</font></h3>
 
<p>
One important special case for a pressure drop is the friction at the
wall of a pipe under the assumption of quasi steady state flow (i.e., the
mass flow rate varies only slowly). In this section it is explained how this case is
handeled in the Modelica_Fluid library for commercial pipes with
<b>nonuniform roughness</b>, including the smooth pipe
as a special case (see 
<a href=\"Modelica://Modelica_Fluid.PressureLosses.WallFrictionAndGravity\">PressureLosses.WallFrictionAndGravity</a>
and 
<a href=\"Modelica://Modelica_Fluid.PressureLosses.Utilities.WallFriction\">PressureLosses.Utilities.WallFriction</a>).
The treatment is non-standard in order to get a 
numerically well-posed description.
</p>
 
<p>
For pipes with circular cross section the pressure drop is computed as:
</p>
 
<pre>
   dp = &lambda;(Re,<font face=\"Symbol\">D</font>)*(L/D)*&rho;*v*|v|/2
      = &lambda;(Re,<font face=\"Symbol\">D</font>)*8*L/(&pi;^2*D^5*&rho;)*m_flow*|m_flow|
      = &lambda;2(Re,<font face=\"Symbol\">D</font>)*k2*sign(m_flow);
 
with
   Re     = |v|*D*&rho;/&eta;
          = |m_flow|*4/(&pi;*D*&eta;)     
   m_flow = A*v*&rho;
   A      = &pi;*(D/2)^2
   &lambda;2     = &lambda;*Re^2
   k2     = L*&eta;^2/(2*D^3*&rho;)
</pre>
 
<p>
where
</p>
<ul>
<li> L is the length of the pipe.</li>
<li> D is the diameter of the pipe. If the pipe has not a 
     circular cross section, D = 4*A/P, where A is the cross section
     area and P is the wetted perimeter.</li>
<li> &lambda; = &lambda;(Re,<font face=\"Symbol\">D</font>) is the \"usual\" wall friction coefficient.</li>
<li> &lambda;2 = &lambda;*Re^2 is the used friction coefficient to get a numerically
     well-posed formulation.</li>
<li> Re = |v|*D*&rho;/&eta; is the Reynolds number.</li>
<li> <font face=\"Symbol\">D</font> = <font face=\"Symbol\">d</font>/D is the relative roughness where
     \"<font face=\"Symbol\">d</font>\" is
     the absolute \"roughness\", i.e., the averaged height of asperities in the pipe
     (<font face=\"Symbol\">d</font> may change over time due to growth of surface asperities during
      service, see <i>[Idelchick 1994, p. 85, Tables 2-1, 2-2])</i>.</li>
<li> &rho; is the upstream density.</li>
<li> &eta; is the upstream dynamic viscosity.</li>
<li> v is the mean velocity.</li>
</ul>
<p>
The first form with &lambda; is used and presented in textbooks,
see \"blue\" curve in the next figure:
</p>
 
<IMG SRC=\"../Images/Components/PipeFriction1.png\" ALT=\"PipeFriction1\">
 
<p>
This form is not suited for a simulation program since 
&lambda; = 64/Re if Re &lt; 2000, i.e., a division by zero occurs for 
zero mass flow rate because Re = 0 in this case. 
More useful for a simulation model is the friction coefficient 
&lambda;2 = &lambda;*Re^2, because &lambda;2 = 64*Re if Re &lt; 2000 and 
therefore no problems for zero mass flow rate occur. 
The characteristic of &lambda;2 is shown in the next figure and is
used in Modelica_Fluid:
</p>
 
<IMG SRC=\"../Images/Components/PipeFriction2.png\" ALT=\"PipeFriction2\">
 
<p>
The pressure loss characteristic is divided into three regions:
</p>
 
<ul>
<li> <b>Region 1</b>: 
     For <b>Re &le; 2000</b>, the flow is <b>laminar</b> and the exact solution of the
     3-dim. Navier-Stokes equations (momentum and mass balance) is used under the
     assumptions of steady flow, constant pressure gradient and constant
     density and viscosity (= Hagen-Poiseuille flow) leading to &lambda;2 = 64*Re.
     Therefore:
     <pre> 
        dp = 128*&eta;*L/(&pi;*D^4*&rho;)*m_flow
     </pre>
</li> 
 
<li> <b>Region 3</b>:
     For <b>Re &ge; 4000</b>, the flow is <b>turbulent</b>.
     Depending on the calculation direction (see \"inverse formulation\"
     below) either of two explicite equations are used. If the pressure drop dp 
     is assumed to be known, &lambda;2 = |dp|/k2. The
     Colebrook-White equation
     <i>[Colebrook 1939; Idelchik 1994, p. 83, eq. (2-9)]</i>:
     <pre>
 
        1/sqrt(&lambda;) = -2*lg( 2.51/(Re*sqrt(&lambda;)) + 0.27*<font face=\"Symbol\">D</font>)
     </pre>
     gives an implicit relationship between Re and &lambda;. 
     Inserting &lambda;2 = &lambda;*Re^2 allows to solve this equation analytically 
     for Re:
     <pre>
 
         Re = -2*sqrt(&lambda;2)*lg(2.51/sqrt(&lambda;2) + 0.27*<font face=\"Symbol\">D</font>)
     </pre>
     Finally, the mass flow rate m_flow is computed from Re via
     m_flow = Re*&pi;*D*&eta;/4*sign(dp).
     These are the <b>red</b> curves in the diagrams above.<br>
     If the mass flow rate is assumed known (and therefore implicitly
     also the Reynolds number), then &lambda;2 is computed by an 
     approximation of the inverse of the Colebrook-White equation 
     <i>[Swamee and Jain 1976;
     Miller 1990, p. 191, eq.(8.4)]</i> adapted to &lambda;2:
     <pre>
 
        &lambda;2 = 0.25*(Re/lg(<font face=\"Symbol\">D</font>/3.7 + 5.74/Re^0.9))^2
     </pre>
     The pressure drop is then computed as dp = k2*&lambda;2*sign(m_flow).
     These are the <b>blue</b> curves in the diagrams above.<br>&nbsp;</li>
 
<li> <b>Region 2</b>:
     For <b>2000 &le; Re &le; 4000</b> there is a transition region between laminar
     and turbulent flow. The value of &lambda;2 depends on more factors as just
     the Reynolds number and the relative roughness, therefore only crude
     approximations are possible in this area.<br>
     The deviation from the laminar region depends on the
     relative roughness. A laminar flow at Re=2000 is only reached for smooth pipes.
     The deviation Reynolds number Re1 is computed according to
     <i>[Samoilenko 1968; Idelchik 1994, p. 81, sect. 2.1.21]</i> as: 
     <pre>
 
        Re1 = 745*e^(if <font face=\"Symbol\">D</font> &le; 0.0065 then 1 else 0.0065/<font face=\"Symbol\">D</font>)
     </pre>
     These are the <b>blue</b> curves in the diagrams above.<br>
     Between Re1=Re1(<font face=\"Symbol\">d</font>/D) and Re2=4000, 
     &lambda;2 is approximated by a cubic
     polynomial in the \"lg(&lambda;2) - lg(Re)\" chart (see figures above) such that the
     first derivative is continuous at these two points. In order to avoid
     the solution of non-linear equations, two different cubic polynomials are used
     for the direct and the inverse formulation. This leads to some discrepancies
     in &lambda;2 (= differences between the red and the blue curves).
     This is acceptable, because the transition region is anyway not
     precisely known since the actual friction coefficient depends on
     additional factors and since the operating points are usually
     not in this region.</li>
</ul>
<p>
The absolute roughness <font face=\"Symbol\">d</font> has usually to
be estimated. In <i>[Idelchik 1994, pp. 105-109,
Table 2-5; Miller 1990, p. 190, Table 8-1]</i> many examples are given.
As a short summary:
</p>
<table border=1 cellspacing=0 cellpadding=2>
  <tr><td><b>Smooth pipes</b></td>
      <td>Drawn brass, coper, aluminium, glass, etc.</td>
      <td><font face=\"Symbol\">d</font> = 0.0025 mm</td>
  </tr>
  <tr><td rowspan=\"3\"><b>Steel pipes</b></td>
      <td>New smooth pipes</td>
      <td><font face=\"Symbol\">d</font> = 0.025 mm</td>
  </tr>
  <tr><td>Mortar lined, average finish</td>
      <td><font face=\"Symbol\">d</font> = 0.1 mm</td>
  </tr>
  <tr><td>Heavy rust</td>
      <td><font face=\"Symbol\">d</font> = 1 mm</td>
  </tr>
  <tr><td rowspan=\"3\"><b>Concrete pipes</b></td>
      <td>Steel forms, first class workmanship</td>
      <td><font face=\"Symbol\">d</font> = 0.025 mm</td>
  </tr>
  <tr><td>Steel forms, average workmanship</td>
      <td><font face=\"Symbol\">d</font> = 0.1 mm</td>
  </tr>
  <tr><td>Block linings</td>
      <td><font face=\"Symbol\">d</font> = 1 mm</td>
  </tr>
</table>
<p>
The equations above are valid for incompressible flow.
They can also be applied for <b>compressible</b> flow up to about <b>Ma = 0.6</b>
(Ma is the Mach number) with a maximum error in &lambda; of about 3 %.
The effect of gas compressibility in a wide region can be taken into
account by the following formula derived by Voronin
<i>[Voronin 1959; Idelchick 1994, p. 97, sect. 2.1.81]</i>:
</p>
<pre>
  &lambda;_comp = &lambda;*(1 + (&kappa;-1)/2 * Ma^2)^(-0.47)
</pre>
<p>
where &kappa; is the isentropic coefficient
(for ideal gases, &kappa; is the ratio of specific heat capacities cp/cv). 
An appreciable decrease in the coefficent \"&lambda;_comp\" is observed
only in a narrow transonic region and also at supersonic flow velocities
by about 15% <i>[Idelchick 1994, p. 97, sect. 2.1.81]</i>.
This effect is not yet included in Modelica_Fluid. 
Another restriction is that the pressure drop model is valid 
only for steady state or slowly changing mass flow rate. 
For large fluid acceleration, the pressure drop depends additionally 
on the frequency of the changing mass flow rate.
</p>
 
<h4><font color=\"#008000\">Inverse formulation</font></h4>
 
<p>
In the \"Advanced menu\" it is possible via parameter
\"from_dp\" to define in which form the
pressure drop equation is actually evaluated (<b>default</b> is from_dp = <b>true</b>):
</p>
<pre>
   from_dp = <b>true</b>:   m_flow = f1(dp)
           = <b>false</b>:  dp     = f2(m_flow)
</pre>
<p>
\"from_dp\" can be useful to avoid nonlinear systems of equations
in cases where the inverse pressure loss function is needed.
</p>
<p>
At the 34th Modelica meeting in Vienna it was discussed to introduce
a language element for alternatives, such that the tool can
figure out what alternative to use. If this would be available,
parameter from_dp could be removed and the equations would
be written as:
</p>
<pre>
  alternative
    // m_flow = f1(dp);
  or
    // dp = f2(m_flow);
  end alternative;
</pre>
<p>
The tool has then \"somehow\" to select the better alternative.
Further research is needed to develop appropriate symbolic
transformation algorithms.
</p>
 
 
<h4><font color=\"#008000\">Summary</font></h4>
 
<p>
A detailed pressure drop model for pipe wall friction is
provided in the form m_flow = f1(dp, <font face=\"Symbol\">D</font>) or
dp = f2(m_flow, <font face=\"Symbol\">D</font>).
These functions are continuous and differentiable, 
are provided in an explicit form without solving non-linear equations, 
and do behave well also at small mass flow rates. This pressure drop
model can be used stand-alone in a static momentum balance and in 
a dynamic momentum balance as the friction pressure drop term. 
It is valid for incompressible and compressible flow up to a Mach number of 0.6.
</p>
 
<h4><font color=\"#008000\">References</font></h4>
 
<dl><dt>Colebrook F. (1939):</dt>
    <dd><b>Turbulent flow in pipes with particular reference to the transition
         region between the smooth and rough pipe laws</b>.
         J. Inst. Civ. Eng. no. 4, 14-25.</dd>
    <dt>Idelchik I.E. (1994):</dt>
    <dd><a href=\"http://www.begellhouse.com/books/00c0f05b040d2ec0.html\"><b>Handbook
        of Hydraulic Resistance</b></a>. 3rd edition, Begell House, ISBN
        0-8493-9908-4</dd>
    <dt>Miller D. S. (1990):</dt>
    <dd><b>Internal flow systems</b>.
    2nd edition. Cranfield:BHRA(Information Services).</dd>
    <dt>Samoilenko L.A. (1968):</dt>
    <dd><b>Investigation of the Hydraulic Resistance of Pipelines in the
        Zone of Transition from Laminar into Turbulent Motion</b>.
        Thesis (Cand. of Technical Science), Leningrad.</dd>
    <dt>Swamee P.K. and Jain A.K. (1976):</dt>
    <dd><b>Explicit equations for pipe-flow problems</b>.
         Proc. ASCE, J.Hydraul. Div., 102 (HY5), pp. 657-664.</dd>
    <dt>Voronin F.S. (1959):</dt>
    <dd><b>Effect of contraction on the friction coefficient in a
           turbulent gas flow</b>.
           Inzh. Fiz. Zh., vol. 2, no. 11, pp. 81-85.</dd>
</dl>
 
</html>
"));
  end WallFriction;
    
  class SemiLinearDefinition "SemiLinear definition" 
      
    annotation (Documentation(info="<html>
<h3><font color=\"#008000\" size=5>SemiLinear definition</font></h3>
<p>
In section <a href=\"Modelica:Modelica_Fluid.UsersGuide.ComponentDefinition.UpstreamDiscretization\">Upstream discretization</a>
the usage of the semiLinear() operator is explained to
describe reversing fluid flow. Below, the updated definition
of the semiLinear() operator is shown, as it will be present in the
next release of the Modelica language specification
(the implementation in Dymola follows the specification and 
the recommendation below and does
not generate an event at zero mass flow rate):
</p>
<table border=1 cellspacing=0 cellpadding=2>
 <tr>
  <td width=160> <pre><b>semiLinear</b>(x,positiveSlope,
             negativeSlope)</pre></td>
  <td>Returns \"smooth(0, if x &gt;= 0 then positiveSlope*x else negativeSlope*x)\". 
      The result is of type Real. For non-scalar arguments the function is vectorized.
      <i>[Note, how smooth() is handeled is a quality
         of implementation, e.g., a tool might or might not generate
         an event when x is crossing zero. Usually, it is more efficient
         to not generate an event]</i></td>
 </tr>
</table>
<p>In some situations, equations with the
semiLinear function become underdetermined if the first argument (x) becomes
zero, i.e., there are an infinite number of solutions. It is <i>recommended</i>
that the following rules are used to transform the equations during the
translation phase in order to select one meaningful solution in such cases:</p>
<p><b>Rule 1</b>: The equations</p>
<p>y = semiLinear(x,<b>sa</b>, s1);<br>
y = semiLinear(x, s1, s2);<br>
y = semiLinear(x, s2, s3);<br>
   ...<br>
y = semiLinear(x, s<sub>N</sub>, <b>sb</b>);<br>
   ....</p>
<p>may be replaced by</p>
<p>s1 = <b>noEvent</b>(<b>if</b> x &gt;= 0 <b>then </b>sa <b>else </b>sb);<br>
s2 = s1;<br>
s3 = s2;<br>
   ...<br>
s<sub>N</sub> = s<sub>N-1</sub>;<br>
y =<b> smooth</b>(0,<b> noEvent</b>(<b>if</b> x &gt;= 0 <b>then</b> sa*x <b>else</b>
sb*x));</p>
<p><br>
Remarks:</p>
<p>
Due to this transformation, the ambiguity for x
== 0 is removed. In the original formulation, there is an infinite number of
solutions for s1, s2,.... In the transformed equation, s1 = sa for x == 0. If
the relation would be changed to \"if x&gt;0\", then, s1 = sb. Therefore, the
ambiguity is resolved by this transformation rule to set s1 = s2 = .... = sa for
x == 0.<br>
A tool might also pick the other solution, i.e., s1 = sb for x == 0. This is
useful if, e.g., x(max=0) is defined.</p>
<p>
The if-expression is discontinuous and therefore
usually an event should be generated. In this particular case, the original
equations have \"si*x\" and \"sj*x\" in the equation. This means that the product
is continuous at x=0 and therefore the noEvent(..) operator can be used here.</p>
<p>&nbsp;</p>
<p><b>Rule 2</b>: The equations</p>
<p>x = 0;<br>
y = 0;<br>
y = semiLinear(x, sa, sb);</p>
<p>may be replaced by</p>
<p>x = 0<br>
y = 0;<br>
sa = sb;</p>
<p>&nbsp;</p>
<p><b>Rule 3</b>: The expression </p>
<p>         semiLinear(x, sa, sb)</p>
<p>may be replaced by x*sa if&nbsp; x.min
&gt;= 0, or by x*sb if x.max &lt;= 0<i> [see comment below]<br>
<br>
</i></p>
<p><i>[For symbolic transformations, the
following property is useful (this follows from the definition):</i></p>
<p>semiLinear(m_flow, port_h, h);</i></p>
<p><i>is identical to </i></p>
<p>-semiLinear(-m_flow, h, port_h);</i></p>
<p><i>Note, in combination with rule 1, the ambiguity
at \"m_flow == 0\" is resolved differently leading to a different result at m_flow
= 0.</i></p>
<p><i>The semiLinear function is designed to
handle reversing flow in fluid systems, such as</i></p>
<p><i>H_flow =semiLinear(m_flow, port.h, h);</i></p>
<p><i>i.e., the enthalpy flow rate H _flow is
computed from the mass flow rate m_flow and the upstream specific enthalpy
depending on the flow direction. In some applications, reversing flow cannot
occur. In such a case, appropriate min/max attributes might be set for m_flow
in connectors, in order to remove the unnecessary if-branch by using Rule 3.
This improves usually both the efficiency and the robustness (provided a
reversal flow does really not occur).</i></p>
<p><i>If 3 pipes are connected together and heat
conduction is neglected, a system of equations of the following form occurs:</i></p>
<pre>      H_flow1 = semiLinear(m_flow1, port.h, h1);
      H_flow2 = semiLinear(m_flow2, port.h, h2);
      H_flow3 = semiLinear(m_flow3, port.h, h3);
            0 = m_flow1 + m_flow2 + m_flow3;
            0 = H_flow1 + H_flow2 + H_flow3;
</pre>
<p><i>where 5 variables of the 10 unknowns are
computed somewhere else, depending on the components that are connected
together. If m_flow1 = m_flow2 = m_flow3 = 0, then it follows from the first
three equations that H_flow1 = H_flow2 = H_flow3 = 0 and therefore, there is
no longer an equation in which port.h appears. Since port.h is also not present
in any other equation of the system, the equations are fulfilled for an
arbitrary value of port.h. In this case (and also for 4 and more pipes that are
connected together), it seems no longer possible to transform the equations to
remove this ambiguity. Instead, during integration the simulator has to detect
this situation and has to compute a meaningful value, e.g., value of port.h in
the last accepted step. Note, If only 2 pipes are connected together, the
ambiguity is removed by rule 2.</i></p>
</html>
"));
  end SemiLinearDefinition;
  end ComponentDefinition;
  
  class ReleaseNotes "Release notes" 
    
    annotation (Documentation(info="<HTML>
<h3><font color=\"#008000\" size=5>Release notes</font></h3>
 
<h3><font color=\"#008000\">Version 1.0 Beta 3, 2007-06-05</font></h3>
 
<p>
Changes according to the Modelica Design Meetings since the
Modelica'2006 conference, especially, improved initialization,
changed Source components (input connectors must be enabled),
improved tank component, moved test models from Examples to
new package Test, many more test models, etc.
This version is slightly non-backward compatible to version 1.0 Beta 2.
</p>


<h3><font color=\"#008000\">Version 1.0 Beta 2, 2006-08-28</font></h3>
 
<p>
Package considerably restructured and some new components added.
New examples (ControlledTankSystem, AST_BatchPlant).
</p>


<h3><font color=\"#008000\">Version 0.96, 2006-01-08</font></h3>
 
<ul>
<li> New package Modelica_Fluid.PressureLosses.</li>
<li> New package Modelica_Fluid.WorkInProgress.</li>
<li> New components in Modelica_Fluid.Components:<br>
     ShortPipe, OpenTank, ValveDiscrete, StaticHead.</li>
<li> New components in Modelica_Fluid.Examples.</li>
<li> Improved users guide.</li>
</ul>
 
<h3><font color=\"#008000\">Version 0.910, 2005-10-25</font></h3>
<ul>
<li> Changes as decided on 41th-45th Modelica Design Meetings
     (details, see minutes).
</ul>
<h3><font color=\"#008000\">Version 0.900, 2004-10-18</font></h3>
<ul>
<li> Changes as decided on 40th Modelica Design Meeting in Dresden 
     (see also minutes)
</ul>
<h3><font color=\"#008000\">Version 0.794, 2004-05-31</font></h3>
<ul>
<li> Sensors.mo, Examples/DrumBoiler.mo: extend sensors with user choice
     for measurement unit.</li>
<li> Components.mo, Types.mo: moved components and types to 
     package Examples.</li>
<li> Moved Examples from <b>file</b> Modelica_Fluid/package.mo to 
     Modelica.Media/Examples <b>subdirectory</b> and created separate 
     file per sub-package. This shall simplify the maintenance of
     examples by different authors</li>
<li> Moved Interfaces from file Modelica_Fluid/package.mo to 
     Modelica_Fluid/Interfaces.mo</li>
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
<li><i>Oct., 2003</i><br>
       by Martin Otter: Adapted to latest design of the Modelica.Media
       library.<br>
       by Ruediger Franke: Included sensor components and 
       Modelica_Fluid.Examples.DrumBoiler example.</li>
<li><i>Sept., 2003</i><br>
       by Martin Otter: Changes according to the decisions of the
       Modelica design meeting in Dearborn, Sept. 2-4, 2003.
       Fluid library splitt in to two packages: Modelica.Media
       that contains the media models and Modelica_Fluid that
       contains fluid flow components. Modelica.Media is
       independent of Modelica_Fluid and my be used also from
       other packages that may have a different design as
       Modelica_Fluid.
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
<li><i>Feb., 2003</i><br>
       by Martin Otter: Included several elementary components and
       a model for moisted air. Some elementary components, such as
       FixedAmbient, are adapted versions from the SimpleFlow fluid library
       of Anton Haumer.</li>
<li><i>Dec., 2002</i><br>
       by Hubertus Tummescheit:
       Improved version of the high precision water model
       (Copy from ThermoFluid library, code reorganization,
       enhanced documentation, additional functions).</li>
<li><i>Nov. 30, 2002</i><br>
       by Martin Otter: Improved the design from the design meeting:
       Adapted to Modelica standard library 1.5,
       added \"choicesAllMatching=true\" annotation,
       added short documentation to \"Interfaces\",
       added packages \"Examples\" and \"Media\" (previously called \"Properties\")
       from previous versions and adapted them to the updated
       \"Interfaces\" package.</li>
<li><i>Nov. 20-21, 2002</i><br>
       by Hilding Elmqvist, Mike Tiller, Allan Watson, John Batteh, Chuck Newman,
       Jonas Eborn: Improved at the 32nd Modelica Design Meeting.</li>
<li><i>Nov. 11, 2002</i><br>
       by Hilding Elmqvist, Martin Otter: improved version.</li>
<li><i>Nov. 6, 2002</i><br>
       by Hilding Elmqvist: first version of the basic design.</li>
</ul>
</HTML>
"));
  equation 
    
  end ReleaseNotes;
  
class Contact "Contact" 
    
    annotation (Documentation(info="<html>
<h3><font color=\"#008000\" size=5>Contact</font></h3>
 
<dl>
<dt><b>Development leader:</b>
<dd>Francesco Casella<br>
    Dipartimento di Elettronica e Informazione<br>
    Politecnico di Milano<br>
    Via Ponzio 34/5<br>
    I-20133 Milano, Italy<br>
    email: <A HREF=\"mailto:casella@elet.polimi.it\">casella@elet.polimi.it</A><br>
</dl>
<p><b>Acknowledgements:</b></p>
<p>
The development of this library has been a collaborative effort 
and many have contributed.
</p>
<ul>
<li> The essential design of the Fluid library, especially
     the component interfaces (Modelica_Fluid.Interfaces.FluidPort)
     to allow arbitrary connections between components  
     that fulfill the balance equations, as well as handling
     reversing flow with the semiLinear() operator,
     is from Hilding Elmqvist.</li>
<li> The Fluid library development was organized in 2002-2004 by Martin 
     Otter and since 2004 it is organized by Francesco Casella.</li>
<li> Francesco Casella included several components of his ThermoPower
     library with some rewriting.</li>
<li> The following people contributed to the fluid component models, 
     examples, and the further design of the library
     (alphabetical list):<br>
     John Batteh,
     Francesco Casella, Jonas Eborn, Hilding Elmqvist, 
     R&uuml;diger Franke, Henning Knigge,
     Sven Erik Mattsson, Chuck Newman, Hans Olsson,
     Martin Otter, Katrin Pr&ouml;l&szlig;,
     Christoph Richter, Mike Tiller, Hubertus Tummescheit,
     Allan Watson.</li>
</ul>
</html>"));
end Contact;
end UsersGuide;
end Modelica_Fluid;
