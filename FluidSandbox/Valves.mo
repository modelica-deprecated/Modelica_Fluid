package Valves "Valve models" 
  extends Icons.VariantLibrary;
  model ValveLinear "Valve for water/steam flows with linear pressure drop" 
    
    extends Interfaces.PartialComponent;
    extends FluidInterface.PartialTransportIsenthalpic;
    
    parameter PowerFluid.Types.HydraulicConductance Kv 
      "Hydraulic conductance at full opening";
    Modelica.Blocks.Interfaces.RealInput opening 
    annotation (extent=[-20,70; 20,110],   rotation=-90);
  equation 
    port_a.m_flow = Kv*opening*(port_a.p - port_b.p);
    
  annotation (
    Icon(
        Polygon(points=[-100,50; -100,-50; 0,0; -100,50],  style(
            color=0,
            thickness=2,
            fillPattern=1)),
        Line(points=[0,60; 0,0],   style(
            color=0,
            thickness=2,
            fillPattern=1)),
        Rectangle(extent=[-20,70; 20,50],   style(
            color=0,
            fillColor=0,
            fillPattern=1)),
        Polygon(points=[100,50; 0,0; 100,-50; 100,50],  style(
            color=0,
            thickness=2,
            fillPattern=1)),
           Text(extent=[-143,-66; 148,-106],  string="%name")));
  end ValveLinear;

  model ValveLinearAA 
    "Valve for water/steam flows with linear pressure drop with two PortA's" 
    
    extends Interfaces.PartialComponent;
    extends FluidInterface.PartialTransportIsenthalpicAA;
    
    parameter PowerFluid.Types.HydraulicConductance Kv 
      "Hydraulic conductance at full opening";
    Modelica.Blocks.Interfaces.RealInput opening 
    annotation (extent=[-20,70; 20,110],   rotation=-90);
  equation 
    port_a.m_flow = Kv*opening*(port_a.p - port_b.p);
    
  annotation (
    Icon(
        Polygon(points=[-100,50; -100,-50; 0,0; -100,50],  style(
            color=0,
            thickness=2,
            fillPattern=1)),
        Line(points=[0,60; 0,0],   style(
            color=0,
            thickness=2,
            fillPattern=1)),
        Rectangle(extent=[-20,70; 20,50],   style(
            color=0,
            fillColor=0,
            fillPattern=1)),
        Polygon(points=[100,50; 0,0; 100,-50; 100,50],  style(
            color=0,
            thickness=2,
            fillPattern=1)),
           Text(extent=[-143,-66; 148,-106],  string="%name")));
  end ValveLinearAA;

  model ValveLinearAB 
    "Valve for water/steam flows with linear pressure drop with a PortA and PortB each" 
    
    extends Interfaces.PartialComponent;
    extends FluidInterface.PartialTransportIsenthalpicAB;
    
    parameter PowerFluid.Types.HydraulicConductance Kv 
      "Hydraulic conductance at full opening";
    Modelica.Blocks.Interfaces.RealInput opening 
    annotation (extent=[-20,70; 20,110],   rotation=-90);
  equation 
    port_a.m_flow = Kv*opening*(port_a.p - port_b.p);
    
  annotation (
    Icon(
        Polygon(points=[-100,50; -100,-50; 0,0; -100,50],  style(
            color=0,
            thickness=2,
            fillPattern=1)),
        Line(points=[0,60; 0,0],   style(
            color=0,
            thickness=2,
            fillPattern=1)),
        Rectangle(extent=[-20,70; 20,50],   style(
            color=0,
            fillColor=0,
            fillPattern=1)),
        Polygon(points=[100,50; 0,0; 100,-50; 100,50],  style(
            color=0,
            thickness=2,
            fillPattern=1)),
           Text(extent=[-143,-66; 148,-106],  string="%name")));
  end ValveLinearAB;
end Valves;
