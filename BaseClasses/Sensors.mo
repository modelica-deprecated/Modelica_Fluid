package Sensors 
  partial model PartialAbsoluteSensor 
    "Partial component to model a sensor that measures a potential variable" 
    
    replaceable package Medium = PackageMedium extends 
      Modelica.Media.Interfaces.PartialMedium "Medium in the sensor" annotation (
        choicesAllMatching =                                                                        true);
    FluidPort_a port(redeclare package Medium = Medium) 
      annotation (extent=[-10,-110; 10,-90],    rotation=90);
    
    annotation (Documentation(info="<html>
<p>
Partial component to model an <b>absolute sensor</b>. Can be used for pressure sensor models.
Use for other properties such as temperature or density is discouraged, because the enthalpy at the connector can have different meanings, depending on the connection topology. Use <tt>PartialFlowSensor</tt> instead.
as signal.
</p>
</html>"),
      Diagram,
      Coordsys(grid=[1,1], scale=0));
  equation 
    port.m_flow = 0;
    port.H_flow = 0;
    port.mXi_flow = zeros(Medium.nXi);
  end PartialAbsoluteSensor;
  
  partial model PartialFlowSensor 
    "Partial component to model sensors that measure flow properties" 
    
    import Modelica.Constants;
    
    replaceable package Medium = PackageMedium extends 
      Modelica.Media.Interfaces.PartialMedium "Medium in the sensor"  annotation (
        choicesAllMatching = true);
    Medium.SpecificEnthalpy h "enthalpy in flow";
    Medium.MassFraction[Medium.nXi] Xi "flow composition";
    
    FluidPort_a port_a(redeclare package Medium = Medium,
                       m_flow(min=if allowFlowReversal then -Constants.inf else 0)) 
      annotation (extent=[-110,-10; -90,10]);
    FluidPort_b port_b(redeclare package Medium = Medium,
                       m_flow(max=if allowFlowReversal then +Constants.inf else 0)) 
      annotation (extent=[110,-10; 90,10]);
    
    parameter Modelica_Fluid.Types.FlowDirectionWithGlobalDefault.Temp 
      flowDirection=
              Modelica_Fluid.Types.FlowDirectionWithGlobalDefault.UseGlobalFluidOption 
      "Unidirectional (port_a -> port_b) or bidirectional flow component" 
       annotation(Dialog(tab="Advanced"));
    
    annotation (Documentation(info="<html>
<p>
Partial component to model a <b>sensor</b> that measures any intensive properties
of a flow, e.g., to get temperature or density in the flow
between fluid connectors.<br>
The model includes zero-volume balance equations. Sensor models inheriting from
this partial class should add a medium instance to calculate the measured property.
</p>
</html>"),
      Diagram,
      Coordsys(grid=[1,1], scale=0));
  protected 
    outer Modelica_Fluid.Components.FluidOptions fluidOptions 
      "Global default options";
    parameter Boolean allowFlowReversal=
       flowDirection == Modelica_Fluid.Types.FlowDirectionWithGlobalDefault.Bidirectional
       or flowDirection == Modelica_Fluid.Types.FlowDirectionWithGlobalDefault.UseGlobalFluidOption
       and fluidOptions.default_flowDirection ==Modelica_Fluid.Types.FlowDirection.Bidirectional 
      "= false, if flow only from port_a to port_b, otherwise reversing flow allowed"
       annotation(Evaluate=true, Hide=true);
  equation 
    port_a.p   = port_b.p;
    // Local *zero-volume* enthalpy and composition
    port_a.H_flow = semiLinear(port_a.m_flow,port_a.h,h);
    port_b.H_flow = semiLinear(port_b.m_flow,port_b.h,h);
    port_a.mXi_flow = semiLinear(port_a.m_flow,port_a.Xi,Xi);
    port_b.mXi_flow = semiLinear(port_b.m_flow,port_b.Xi,Xi);
    // Static balances
    0 = port_a.m_flow + port_b.m_flow;
    0 = port_a.H_flow + port_b.H_flow;
    zeros(Medium.nXi) = port_a.mXi_flow + port_b.mXi_flow;
  end PartialFlowSensor;
  
protected 
  partial model PartialRelativeSensor 
    "Partial component to model a sensor that measures the difference of effort variables at two ports" 
    
    replaceable package Medium = PackageMedium extends 
      Modelica.Media.Interfaces.PartialMedium "Medium in the sensor"  annotation (
        choicesAllMatching =                                                                         true);
    
    FluidPort_a port_a(redeclare package Medium = Medium) 
      annotation (extent=[-120, -10; -100, 10]);
    FluidPort_b port_b(redeclare package Medium = Medium) 
      annotation (extent=[120, -10; 100, 10]);
    
    annotation (Documentation(info="<html>
<p>
Partial component to model a <b>sensor</b> that measures
the <b>difference between two effort variables</b>, e.g. to obtain the temperature difference 
between fluid connectors.
</p>
</html>"));
  equation 
    port_a.m_flow = 0;
    port_a.H_flow = 0;
    port_a.mXi_flow = zeros(Medium.nXi);
    
    port_b.m_flow = 0;
    port_b.H_flow = 0;
    port_b.mXi_flow = zeros(Medium.nXi);
  end PartialRelativeSensor;
  
end Sensors;
