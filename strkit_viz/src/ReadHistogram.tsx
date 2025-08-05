import { useEffect, useRef } from "react";
import * as d3 from "d3";
import { type Bin } from "d3";
import type { Peaks } from "./types";

const palettes = [d3.interpolateBlues, d3.interpolateReds];

const normalPdf = (mu: number, sigma: number) => (x: number) =>
  (1 / (sigma * Math.sqrt(2 * Math.PI))) * Math.pow(Math.E, -0.5 * Math.pow((x - mu) / sigma, 2));

function Histogram(data: number[], {
  // normal stuff:
  peaks,

  // hist stuff:
  x, // given d in data, returns the (quantitative) x-value
  y = () => 1, // given d in data, returns the (quantitative) weight
  thresholds = 40, // approximate number of bins to generate, or threshold function
  marginTop = 20, // top margin, in pixels
  marginRight = 30, // right margin, in pixels
  marginBottom = 30, // bottom margin, in pixels
  marginLeft = 40, // left margin, in pixels
  width = 640, // outer width of chart, in pixels
  height = 400, // outer height of chart, in pixels
  insetLeft = 0.5, // inset left edge of bar
  insetRight = 0.5, // inset right edge of bar
  xRange = [marginLeft, width - marginRight], // [left, right]
  xLabel, // a label for the x-axis
  yRange = [height - marginBottom, marginTop], // [bottom, top]
  yLabel = "↑ Frequency", // a label for the y-axis
  color = "currentColor", // bar fill color
  showGaussians = true,
}: {
  peaks: Peaks,
  x: (d: number) => number,
  y?: (d: number) => number,
  thresholds: number,
  // --------------------------
  marginTop?: number,
  marginRight?: number,
  marginBottom?: number,
  marginLeft?: number,
  width?: number,
  height?: number,
  insetLeft?: number,
  insetRight?: number,
  // --------------------------
  xLabel: string,
  xRange?: [number, number],
  yRange?: [number, number],
  yLabel: string,
  color: string,
  showGaussians: boolean,
}) {
  // Compute values.
  const X = d3.map(data, x);
  const Y0 = d3.map(data, y);
  const I = d3.range(X.length);

  // Compute bins.
  const bins = d3.bin().thresholds(thresholds).value(i => X[i])(I);
  const Y = Array.from(bins, I => d3.sum(I, i => Y0[i]));

  // Compute default domains.
  const xDomain = [bins[0].x0 - 1, bins[bins.length - 1].x1 + 1];
  const yDomain = [0, d3.max(Y)];

  // Construct scales and axes.
  const xScale = d3.scaleLinear(xDomain, xRange);
  const yScale = d3.scaleLinear(yDomain, yRange);
  const xAxis = d3.axisBottom(xScale).ticks(Math.min(width / 80, bins.length + 1)).tickSizeOuter(0);
  const yAxis = d3.axisLeft(yScale).ticks(height / 40);
  const yFormat = yScale.tickFormat(100);

  const svg = d3.create("svg")
    .attr("width", width)
    .attr("height", height)
    .attr("viewBox", [0, 0, width, height])
    .attr("style", "max-width: 100%; height: auto; height: intrinsic;");

  svg.append("g")
    .attr("transform", `translate(${marginLeft},0)`)
    .call(yAxis)
    .call(g => g.select(".domain").remove())
    .call(g => g.selectAll(".tick line").clone()
      .attr("x2", width - marginLeft - marginRight)
      .attr("stroke-opacity", 0.1))
    .call(g => g.append("text")
      .attr("x", -marginLeft)
      .attr("y", 10)
      .attr("fill", "currentColor")
      .attr("text-anchor", "start")
      .text(yLabel));

  const barWidth = (d: Bin<number, number>) => Math.max(0, xScale(d.x1) - xScale(d.x0) - insetLeft - insetRight);

  svg.append("g")
    .attr("fill", color)
    .selectAll("rect")
    .data(bins)
    .join("rect")
    .attr("x", d => xScale(d.x0) + insetLeft - barWidth(d)/2)
    .attr("width", barWidth)
    .attr("y", (d, i) => yScale(Y[i]))
    .attr("height", (d, i) => yScale(0) - yScale(Y[i]))
    .append("title")
    .text((d, i) => [`${d.x0} ≤ x < ${d.x1}`, yFormat(Y[i])].join("\n"));

  if (peaks && showGaussians) {
    const xs = [];
    for (let i = xDomain[0]; i < xDomain[1]; i += 0.01) {
      xs.push(i);
    }
    const line = d3.line()
      .curve(d3.curveBasis)
      .x(d => xScale(d[0]))
      .y(d => yScale(d[1]));
    const max = yDomain[1];
    const pdfMaxes = [];
    for (let i = 0; i < peaks.modal_n; i++) {
      const stDev = peaks.stdevs[i] || 0.05;
      pdfMaxes.push(1 / (Math.sqrt(2 * Math.PI) * stDev));
    }
    const maxPdfMax = Math.max(...pdfMaxes);
    for (let i = 0; i < peaks.modal_n; i++) {
      console.log(peaks.stdevs);
      const stDev = peaks.stdevs[i] || 0.05;
      const pdf = normalPdf(peaks.means[i], stDev);
      svg.append("path")
        .datum(xs.map(xx => [xx, pdf(xx) * (max / maxPdfMax)])) // peaks.weights[i]
        .attr("fill", "none")
        .attr("stroke", palettes[i](0.7))
        .attr("stroke-width", 2.5)
        .attr("stroke-linejoin", "round")
        // @ts-ignore
        .attr("d", line);
    }
  }

  svg.append("g")
    .attr("transform", `translate(0,${height - marginBottom})`)
    .call(xAxis)
    .call(g => g.append("text")
      .attr("x", width - marginRight)
      .attr("y", 27)
      .attr("fill", "currentColor")
      .attr("text-anchor", "end")
      .text(xLabel));

  return svg.node();
}

const ReadHistogram = ({ cns, peaks }: { cns: number[], peaks: Peaks }) => {
  const svg = useRef(null);

  console.log(cns, cns.length);

  const thresholds = Math.min(100, Math.max(...cns) - Math.min(...cns));

  const hist = Histogram(cns, {
    peaks,
    thresholds,
    x: d => d,
    xLabel: "Copy number →",
    yLabel: "↑ Read count",
    width: 800,
    height: 500,
    color: "#9FA199",
    showGaussians: true, // ui.showGaussians,
  });

  useEffect(() => {
    if (svg.current) {
      svg.current.appendChild(hist);
    }
  }, []);

  return <div ref={svg} />;
};

export default ReadHistogram;
