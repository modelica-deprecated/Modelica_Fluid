package Modelica_Fluid "Fluid package that should be included into package Modelica"
annotation (
  version="0.900",
  versionDate="2004-10-18",
  preferedView="info",
  Settings(NewStateSelection=true),
  uses( Modelica(version="2.2")),
  Documentation(info="<html>
<p>
This library provides basic components and property models
to model <b>1-dimensional thermo-fluid flow</b> systems of
a <b>single substance</b> or of a <b>mixture of substances</b> 
with optional <b>multiple phases</b>. The goal is to include 
this library in the Modelica standard library. The Modelica_Fluid
library uses the media models from the Modelica.Media library
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
a standardized interface to components describing
1-dimensional fluid flow in networks of pipes. A unique feature is that the
component equations and the media models are decoupled.
All components are implemented such that they can be used for
media from the Modelica.Media library. This means especially that an
incompressible or compressible medium, a single or a multiple
substance medium with one or more phases might be used for 
every component model in the Modelica_Fluid library.
</p>
<p>
This users guide is just a start and will be improved
considerably.
</p>
</HTML>"));
  
  model ReleaseNotes "Release notes" 
    
    annotation (Documentation(info="<HTML>
<h3><font color=\"#008000\" size=5>Release notes</font></h3>
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
       Fluid library splitt in to two packages: Modelica.Media
       that contains the media models and Modelica_Fluid that
       contains fluid flow components. Modelica.Media is
       independent of Modelica_Fluid and my be used also from
       other packages that may have a different design as
       Modelica_Fluid.
<li><i>Oct., 2003</i><br>
       by Martin Otter: Adapted to latest design of the Modelica.Media
       library.<br>
       by Ruediger Franke: Included sensor components and 
       Modelica_Fluid.Examples.DrumBoiler example.</li>
</ul>
</HTML>
"));
  equation 
    
  end ReleaseNotes;
end UsersGuide;

replaceable package PackageMedium = Modelica.Media.Interfaces.PartialMedium 
  "To allow change of default medium for all components" annotation (
    choicesAllMatching=true);
end Modelica_Fluid;
