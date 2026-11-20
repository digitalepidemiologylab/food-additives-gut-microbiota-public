# Food Additive Network Visualization

## Overview
This interactive visualization tool displays food additive co-occurrence networks, showing which additives frequently appear together in food products. The network represents additives as nodes (colored by category) and their co-occurrences as edges.

## Required Files

Two input files are needed to generate the visualization:

| File | Location | Description |
|------|----------|-------------|
| Annotation CSV | `data/additive_name_mapping_v3.csv` | Contains E-code identifiers, additive names, and category assignments |
| Co-occurrence JSON | `data/additive_analysis_files/json/` | Contains additive co-occurrence data for specific food groups |

## Usage

1. Open the HTML file in a web browser.

2. Upload the required files using the corresponding buttons:
   - **"Upload Network JSON File"**: Select a co-occurrence JSON file from `data/additive_analysis_files/json/`
   - **"Upload Annotation CSV File"**: Select `additive_name_mapping_v3.csv` from `data/`

3. Click **"Visualize Network"** to generate the visualization.

## Visualization Features

- **Node size**: Proportional to additive occurrence frequency
- **Node color**: Indicates additive category (see legend)
- **Edge thickness**: Reflects co-occurrence frequency
- **Interactive controls**:
  - Zoom in/out and reset view
  - Toggle node labels on/off
  - Toggle between E-code and additive name display
  - Drag nodes to rearrange layout
  - Hover over nodes for details
  - Click nodes to view co-occurrence information
  - Export visualization as PNG image

## File Formats

### Annotation CSV
```csv
ecode,name,category
e100,Curcumin,colorants
e101,Riboflavin,colorants
```

### Co-occurrence JSON
Array of objects containing additive combinations and their counts within a food group.