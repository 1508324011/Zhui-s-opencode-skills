---
name: string-app
description: STRING database integration and protein-protein interaction network analysis. Use for querying protein interaction networks, enrichment analysis, functional annotation, and network topology analysis. Best for discovering RBP interaction partners, pathway analysis, and functional module identification. Can be used standalone or integrated with Cytoscape.
license: CC-BY
compatibility: opencode
metadata:
  category: protein-interaction-networks
  language: python
  skill-author: K-Dense Inc.
---

# STRING-APP: Protein-Protein Interaction Network Analysis

## Overview

STRING (Search Tool for the Retrieval of Interacting Genes/Proteins) is a comprehensive database of known and predicted protein-protein interactions (PPIs). The string-app skill provides programmatic access to STRING data for network analysis, enrichment studies, and functional annotation. This skill integrates with the string-database skill but focuses on network analysis and Cytoscape integration.

**Current Version**: STRING 12.0
**Coverage**: 24.5 million proteins from 12,000+ organisms
**Interaction Types**: Physical binding, co-expression, co-occurrence, text mining, databases

## When to Use This Skill

Use STRING-APP when you need to:

- **Discover RBP interaction partners** and binding proteins
- **Analyze protein interaction networks** in context of pathways
- **Perform functional enrichment analysis** on gene sets
- **Identify functional modules** within protein networks
- **Predict novel protein functions** based on network neighbors
- **Export networks to Cytoscape** for advanced visualization
- **Query interaction evidence** and confidence scores
- **Compare networks across species**

## Installation and Setup

```bash
# Using pip
pip install stringdb

# Alternative: requests for direct API access
pip install requests pandas networkx
```

## Core Capabilities

### 1. Interaction Types in STRING

**Physical Interactions**:
- Experimentally determined (yeast two-hybrid, co-IP, etc.)
- Protein 3D structures (PDB)

**Functional Associations**:
- Co-expression across experiments
- Co-occurrence in genomes
- Text mining of PubMed abstracts
- Database imports (KEGG, Reactome, etc.)

**Confidence Scoring**: 0-1000 (higher = more confident)

### 2. Network Analysis Features

- **Network expansion**: Find interactors of a protein set
- **Clustering**: Identify functional modules (k-means, MCL)
- **Enrichment**: GO, KEGG, Pfam, InterPro enrichment
- **Topology**: Degree, betweenness, clustering analysis
- **Comparison**: Compare networks between conditions

### 3. Data Export Formats

- **TSV/CSV**: Tabular interaction data
- **Network formats**: Cytoscape, BioPAX, PSI-MI, GML
- **Image**: PNG, SVG network visualizations

## Using This Skill

### Example 1: Query RBP Interaction Network

```python
import requests
import pandas as pd
import networkx as nx

# STRING API endpoint
string_api_url = "https://string-db.org/api"
output_format = "json"
method = "network"

# Your RBP list
rbp_list = ["FUS", "TARDBP", "HNRNPA1", "PABPC1", "ELAVL1", "RBFOX2"]
species = "9606"  # Human NCBI taxonomy ID

# Construct request
request_url = "/".join([string_api_url, output_format, method])
params = {
    "identifiers": "%0d".join(rbp_list),
    "species": species,
    "network_flavor": "evidence",  # confidence, evidence, actions
    "caller_identity": "opencode_user"
}

# Get network data
response = requests.post(request_url, data=params)
network_data = response.json()

# Parse into DataFrame
interactions = []
for edge in network_data:
    interactions.append({
        'protein1': edge['preferredName_A'],
        'protein2': edge['preferredName_B'],
        'score': edge['score'],
        'evidence': edge.get('evidence', '')
    })

df = pd.DataFrame(interactions)
print(f"Found {len(df)} interactions")
print(df.head(10))
```

### Example 2: Create NetworkX Graph

```python
# Create network from STRING data
G = nx.Graph()

# Add edges with confidence scores
for edge in network_data:
    p1 = edge['preferredName_A']
    p2 = edge['preferredName_B']
    score = edge['score'] / 1000  # Normalize to 0-1
    
    G.add_edge(p1, p2, weight=score, evidence=edge.get('evidence', ''))

# Add node attributes
for protein in rbp_list:
    if protein in G:
        G.nodes[protein]['type'] = 'RBP'
        G.nodes[protein]['query_protein'] = True

# Network statistics
print(f"Nodes: {G.number_of_nodes()}")
print(f"Edges: {G.number_of_edges()}")
print(f"Density: {nx.density(G):.3f}")

# Degree analysis
degrees = dict(G.degree())
sorted_degree = sorted(degrees.items(), key=lambda x: x[1], reverse=True)
print("\nTop hub proteins:")
for protein, degree in sorted_degree[:10]:
    print(f"  {protein}: {degree} connections")
```

### Example 3: Functional Enrichment Analysis

```python
# Enrichment analysis
enrichment_method = "enrichment"
enrichment_url = "/".join([string_api_url, output_format, enrichment_method])

enrichment_params = {
    "identifiers": "%0d".join(rbp_list),
    "species": species,
    "caller_identity": "opencode_user"
}

enrichment_response = requests.post(enrichment_url, data=enrichment_params)
enrichment_data = enrichment_response.json()

# Parse enrichment results
enrichment_df = pd.DataFrame(enrichment_data)

# Filter significant results (FDR < 0.05)
significant = enrichment_df[enrichment_df['fdr'] < 0.05].sort_values('fdr')

print("Top enriched pathways:")
print(significant[['category', 'term', 'description', 'fdr']].head(10))

# Category breakdown
categories = enrichment_df['category'].value_counts()
print("\nEnrichment categories:")
print(categories)
```

### Example 4: Network Expansion

```python
# Find first neighbors (interaction partners)
interaction_partners_method = "interaction_partners"
partners_url = "/".join([string_api_url, output_format, interaction_partners_method])

partners_params = {
    "identifiers": "%0d".join(rbp_list),
    "species": species,
    "limit": 50,  # Top 50 interactors per protein
    "caller_identity": "opencode_user"
}

partners_response = requests.post(partners_url, data=partners_params)
partners_data = partners_response.json()

# Build expanded network
G_expanded = nx.Graph()

for node in partners_data:
    query_protein = node['queryItem']
    for partner in node['interactors']:
        partner_name = partner['preferredName']
        score = partner['score'] / 1000
        
        G_expanded.add_edge(query_protein, partner_name, weight=score)
        
        # Mark query proteins
        G_expanded.nodes[query_protein]['query'] = True

print(f"Expanded network: {G_expanded.number_of_nodes()} nodes, {G_expanded.number_of_edges()} edges")

# Identify novel interactors
novel_interactors = set(G_expanded.nodes()) - set(rbp_list)
print(f"\nNovel interacting proteins: {len(novel_interactors)}")
```

### Example 5: Export to Cytoscape

```python
# Export network for Cytoscape
nx.write_graphml(G_expanded, 'rbp_string_network.graphml')

# Create node attributes file
node_attrs = []
for node in G_expanded.nodes():
    is_rbp = node in rbp_list
    degree = G_expanded.degree(node)
    betweenness = nx.betweenness_centrality(G_expanded)[node]
    
    node_attrs.append({
        'node': node,
        'is_query_rbp': is_rbp,
        'degree': degree,
        'betweenness': betweenness
    })

node_df = pd.DataFrame(node_attrs)
node_df.to_csv('rbp_network_node_attributes.csv', index=False)

print("Files exported:")
print("  - rbp_string_network.graphml (network)")
print("  - rbp_network_node_attributes.csv (attributes)")
print("\nOpen in Cytoscape: File → Import → Network from File")
```

### Example 6: Clustering and Module Detection

```python
from sklearn.cluster import KMeans
import numpy as np

# Calculate network proximity matrix
proteins = list(G_expanded.nodes())
n = len(proteins)
proximity_matrix = np.zeros((n, n))

# Use shortest path length as distance
for i, p1 in enumerate(proteins):
    for j, p2 in enumerate(proteins):
        if i != j:
            try:
                dist = nx.shortest_path_length(G_expanded, p1, p2)
                proximity_matrix[i, j] = 1 / (dist + 1)
            except nx.NetworkXNoPath:
                proximity_matrix[i, j] = 0

# K-means clustering
kmeans = KMeans(n_clusters=4, random_state=42)
clusters = kmeans.fit_predict(proximity_matrix)

# Assign clusters
for protein, cluster in zip(proteins, clusters):
    G_expanded.nodes[protein]['cluster'] = int(cluster)

# Analyze clusters
for c in range(4):
    cluster_proteins = [p for p in proteins if G_expanded.nodes[p].get('cluster') == c]
    print(f"\nCluster {c} ({len(cluster_proteins)} proteins):")
    print(f"  {', '.join(cluster_proteins[:10])}")
```

## Integration with Other Skills

- **string-database**: Direct database queries
- **cytoscape**: Advanced network visualization
- **networkx**: Network analysis algorithms
- **pandas**: Data manipulation
- **scikit-learn**: Machine learning for clustering
- **gget**: Gene information queries

## Best Practices

### Network Construction
1. **Confidence threshold**: Filter edges with score < 400
2. **Query size**: Optimal query size is 10-500 proteins
3. **Species selection**: Use correct NCBI taxonomy ID
4. **Evidence types**: Consider interaction evidence types

### Analysis Guidelines
1. **Statistical testing**: Use STRING's built-in enrichment tests
2. **Multiple testing**: Always correct for multiple comparisons
3. **Validation**: Cross-check with experimental data
4. **Interpretation**: Distinguish physical vs. functional interactions

### Visualization
1. **Node sizing**: By degree or centrality
2. **Edge styling**: Thickness by confidence score
3. **Color coding**: By functional categories
4. **Layout**: Force-directed for biological networks

## Advanced Topics

### Custom Evidence Scoring

```python
# Weight different evidence types
def calculate_custom_score(edge_data):
    weights = {
        'experimental': 0.4,
        'database': 0.3,
        'textmining': 0.1,
        'coexpression': 0.1,
        'cooccurrence': 0.1
    }
    
    score = 0
    for evidence, weight in weights.items():
        if evidence in edge_data.get('evidence', ''):
            score += weight
    
    return score
```

### Cross-Species Analysis

```python
# Map proteins between species
species_list = ["9606", "10090", "10116"]  # Human, Mouse, Rat
ortholog_method = "get_orthologs"

for target_species in species_list:
    params = {
        "identifiers": "%0d".join(rbp_list),
        "species": "9606",
        "species_bis": target_species
    }
    # Query orthologs
```

## Troubleshooting

**Issue**: API rate limiting
- **Solution**: Add delays between requests; use bulk queries

**Issue**: Protein ID not found
- **Solution**: Use official gene symbols or Uniprot IDs

**Issue**: Large network exports
- **Solution**: Filter by confidence score before export

## Resources

- **STRING Website**: https://string-db.org/
- **API Documentation**: https://string-db.org/cgi/help?subpage=api
- **Paper**: Szklarczyk et al. (2023) Nucleic Acids Research
- **Cytoscape App**: https://apps.cytoscape.org/apps/stringapp

## Summary

STRING-APP provides powerful capabilities for exploring protein interaction networks, with particular utility for discovering RBP interaction partners and functional context. Its comprehensive database and flexible API make it an essential tool for network-based systems biology research.
