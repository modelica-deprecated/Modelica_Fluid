within ;
package FluidSandbox "Library to study different numerical solution methods in the fluid dynamics domain"
  import SI = Modelica.SIunits;

  annotation (
    version="0.44",
    versionDate="2008-04-03",
    preferedView="info",
    classOrder={"Examples","FluidInterfaces","FluidDiscretizations",
        "*","ConnectionSemantics", "Interfaces","Icons"},
    uses(Modelica_Fluid(version="1.0 Beta 3"),
    Modelica(version="2.2.2")),
    Documentation(info="<html>
<h1>Fluid Sandbox</h1>
<p>
Even though it was proven by various entities that problems of the thermo fluid domain can be modelled successfully using Modelica, it has been somewhat challenging to converge towards a \"one-fits-all\" solution. Consequently, different approaches are frequently discussed. This library is meant as a first step towards a kind of prototyping environment to implement and test connector designs, discretizations and the like. Hopefully, this allows a more rigorous comparison of the different concepts and helps our iterative process.
</p>
<p>
In parallel, an effort was started to summarize some of the options in the fluid domain, and to discuss their respective advantages and drawbacks. 
Turn to <a href=\"https://trac.modelica.org/Modelica_Fluid/wiki/FluidConcepts\">the web page</a> to read and contribute (log in to get access).
</p>
<p>
The FluidConcept library by Katrin Pr&ouml;l&szlig; provided relevant input for this work. 
</p>
<p>
Most component models are based on Modelica_Fluid code.
</p>
</html>", revisions="<html>
<ul>
<li><i>January 2008</i>
   by <a href=\"mailto:michael.sielemann@dlr.de\">Michael Sielemann</a>:<br>
   Created.</li>
</ul>
</html>"));
end FluidSandbox;
