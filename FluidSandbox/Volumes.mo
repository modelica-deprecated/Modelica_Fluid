within FluidSandbox;
package Volumes "Generic mixing volume and other volume type components" 
  extends Icons.VariantLibrary;
  
  model Volume "Volume with inlet and outlet ports (flow reversal is allowed)" 
    
    extends Interfaces.PartialComponent;
    
    extends FluidInterface.PartialLumpedVolume(V_lumped=V, Ws_flow=0);
    
    // Constant volume
    parameter SI.Volume V "Volume";
    
    // Heat port
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a thermalPort 
      "Thermal port" annotation (extent=[-20,88; 20,108]);
    
  equation 
    thermalPort.T = medium.T;
    Qs_flow = thermalPort.Q_flow;
    
    annotation (
      Icon(
        Ellipse(extent=[-100,100; 100,-100], style(
            color=0,
            rgbcolor={0,0,0},
            gradient=3,
            fillColor=68,
            rgbfillColor={170,213,255})),
        Text(extent=[-144,178; 146,116], string="%name"),
        Text(
          extent=[-130,-108; 144,-150],
          style(color=0),
          string="V=%V")),
      Documentation(info="<html>
Ideally mixed volume of constant size with two fluid ports and one medium model. The flow properties are computed from the upstream quantities, pressures are equal in both nodes and the medium model. Heat transfer through a thermal port is possible, it equals zero if the port remains unconnected. The thermal port temperature is equal to the medium temperature.
</html>"),
      Diagram);
  end Volume;
end Volumes;
