package HeatTransfer 
  partial model PartialPipeHeatTransfer 
    replaceable package Medium=Modelica.Media.Interfaces.PartialMedium annotation(Dialog(tab="No input", enable=false));
    parameter Integer n(min=1)=1 "Number of pipe segments" annotation(Dialog(tab="No input", enable=false));
    SI.HeatFlowRate[n] Q_flow "Heat flow rates";
    parameter SI.Area A_h "Total heat transfer area" annotation(Dialog(tab="No input", enable=false));
    parameter SI.Length d_h "Hydraulic diameter" annotation(Dialog(tab="No input", enable=false));
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a[n] thermalPort 
      "Thermal port" 
      annotation (extent=[-20,60; 20,80]);
    SI.Temperature[n] T;
  equation 
    
    annotation (Icon(Ellipse(extent=[-60,64; 60,-56], style(
            color=42,
            rgbcolor={127,0,0},
            gradient=3,
            fillColor=1,
            rgbfillColor={232,0,0})), Text(
          extent=[-38,26; 40,-14],
          style(
            color=42,
            rgbcolor={127,0,0},
            gradient=3,
            fillColor=1,
            rgbfillColor={232,0,0},
            fillPattern=7),
          string="%name")), Documentation(info="<html>
Base class for heat transfer models that can be used in model <b>Pipe</b>.
</html>"));
  end PartialPipeHeatTransfer;
  
  model PipeHT_constAlpha 
    extends PartialPipeHeatTransfer;
    parameter SI.CoefficientOfHeatTransfer alpha0=200;
    annotation(structurallyIncomplete=true, Documentation(info="<html>
Simple heat transfer correlation with constant heat transfer coefficient
</html>"));
  equation 
    for i in 1:n loop
      thermalPort[i].Q_flow=alpha0*A_h/n*(thermalPort[i].T-T[i]);
    end for;
    thermalPort.Q_flow=Q_flow;
  end PipeHT_constAlpha;
  annotation (Documentation(info="<html>
Heat transfer correlations for pipe models
</html>"));
end HeatTransfer;
