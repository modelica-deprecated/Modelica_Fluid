package Sources 
partial model PartialSource "Partial component source with one fluid connector" 
    import Modelica.Constants;
  replaceable package Medium = PackageMedium extends 
      Modelica.Media.Interfaces.PartialMedium "Medium model within the source" 
     annotation (choicesAllMatching=true);
  Interfaces.FluidPort_b port(redeclare package Medium = Medium,
                   m_flow(min=if allowFlowReversal then -Constants.inf else 0)) 
    annotation (extent=[90,-10; 110,10],    rotation=0);
  Medium.BaseProperties medium "Medium in the source";
  parameter Types.FlowDirection.Temp flowDirection=
                   Types.FlowDirection.Unidirectional 
      "Unidirectional (out of port_b) or bidirectional flow component" 
                                                              annotation(Dialog(tab="Advanced"));
  protected 
    parameter Boolean allowFlowReversal=
     flowDirection == Modelica_Fluid.Types.FlowDirection.Bidirectional 
      "= false, if flow only out of port_b, otherwise reversing flow allowed" 
     annotation(Evaluate=true, Hide=true);
equation 
  port.p = medium.p;
  port.H_flow = semiLinear(port.m_flow, port.h, medium.h);
  port.mXi_flow = semiLinear(port.m_flow, port.Xi, medium.Xi);
  annotation (Documentation(info="<html>
<p>
Partial component to model the <b>volume interface</b> of a <b>source</b>
component, such as a mass flow source. The essential
features are:
</p>
<ul>
<li> The pressure in the connection port (= port.p) is identical to the
     pressure in the volume (= medium.p).</li>
<li> The enthalpy flow rate (= port.H_flow) and the mass flow rates of the
     substances (= port.mX_flow) depend on the direction of the mass flow rate.</li>
</ul>
</html>"),
    Diagram,
    Coordsys(grid=[1,1], scale=0));
end PartialSource;
end Sources;
