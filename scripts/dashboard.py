#!/usr/bin/env python
# dashboard.py

import os
import pandas as pd
import streamlit as st
import plotly.express as px
import plotly.graph_objects as go
import numpy as np

st.set_page_config(page_title="COVID-19 Mutation Dashboard", layout="wide")

# Title and description
st.title("COVID-19 Mutation Analysis Dashboard")
st.write("""
This dashboard visualizes mutations detected in COVID-19 samples from different countries
compared to a reference genome sequence.
""")

# Load mutation data
@st.cache_data
def load_data():
    mutations_file = "results/all_mutations.csv"
    if os.path.exists(mutations_file):
        return pd.read_csv(mutations_file)
    else:
        st.error("Mutation data file not found. Please run the analysis pipeline first.")
        return None

data = load_data()

if data is not None:
    # Sidebar filters
    st.sidebar.header("Filters")
    
    countries = ["All"] + sorted(data["Country"].unique().tolist())
    selected_country = st.sidebar.selectbox("Select Country", countries)
    
    mutation_types = ["All"] + sorted(data["Type"].unique().tolist())
    selected_type = st.sidebar.selectbox("Mutation Type", mutation_types)
    
    # Filter data based on selections
    filtered_data = data.copy()
    if selected_country != "All":
        filtered_data = filtered_data[filtered_data["Country"] == selected_country]
    if selected_type != "All":
        filtered_data = filtered_data[filtered_data["Type"] == selected_type]
    
    # Main content in tabs
    tab1, tab2, tab3 = st.tabs(["Overview", "Mutation Analysis", "Country Comparison"])
    
    with tab1:
        col1, col2 = st.columns(2)
        
        with col1:
            st.subheader("Mutation Distribution by Type")
            fig = px.pie(data, names="Type", hole=0.4)
            st.plotly_chart(fig, use_container_width=True)
        
        with col2:
            st.subheader("Mutations by Country")
            country_counts = data.groupby("Country").size().reset_index(name="count")
            fig = px.bar(country_counts, x="Country", y="count")
            st.plotly_chart(fig, use_container_width=True)
        
        # Summary statistics
        st.subheader("Summary Statistics")
        col1, col2, col3, col4 = st.columns(4)
        col1.metric("Total Mutations", len(data))
        col2.metric("Countries", data["Country"].nunique())
        col3.metric("SNPs", len(data[data["Type"] == "SNP"]) if "SNP" in data["Type"].values else 0)
        col4.metric("Structural Variants", 
                   len(data[data["Type"].isin(["INS", "DEL"])]) if any(t in data["Type"].values for t in ["INS", "DEL"]) else 0)
    
    with tab2:
        st.subheader("Mutation Positions Along the Genome")
        
        # Create position bins for better visualization
        max_pos = data["Position"].max()
        bin_size = max(1, max_pos // 100)  # Divide genome into ~100 bins
        data["Position_bin"] = (data["Position"] // bin_size) * bin_size
        
        position_counts = filtered_data.groupby("Position_bin").size().reset_index(name="count")
        fig = px.line(position_counts, x="Position_bin", y="count")
        fig.update_layout(xaxis_title="Genome Position", yaxis_title="Mutation Count")
        st.plotly_chart(fig, use_container_width=True)
        
        # Top mutation positions
        st.subheader("Top Mutation Positions")
        top_positions = filtered_data.groupby("Position").size().reset_index(name="count")
        top_positions = top_positions.sort_values("count", ascending=False).head(10)
        fig = px.bar(top_positions, x="Position", y="count")
        st.plotly_chart(fig, use_container_width=True)
        
        if "SNP" in filtered_data["Type"].values:
            st.subheader("SNP Substitution Patterns")
            snps = filtered_data[filtered_data["Type"] == "SNP"]
            
            # Create a cross-tabulation of reference and alternate bases
            base_counts = pd.crosstab(snps["Ref"], snps["Alt"])
            
            # Fill missing values with 0
            all_bases = sorted(set(snps["Ref"].unique()) | set(snps["Alt"].unique()))
            for base in all_bases:
                if base not in base_counts.index:
                    base_counts.loc[base] = 0
                if base not in base_counts.columns:
                    base_counts[base] = 0
            
            base_counts = base_counts.reindex(sorted(base_counts.index), axis=0)
            base_counts = base_counts.reindex(sorted(base_counts.columns), axis=1)
            
            fig = px.imshow(base_counts, text_auto=True, aspect="auto",
                           labels=dict(x="Alternate Base", y="Reference Base", color="Count"))
            st.plotly_chart(fig, use_container_width=True)
    
    with tab3:
        st.subheader("Mutation Patterns Across Countries")
        
        # Create a mutation signature for each country
        data["Mutation"] = data.apply(lambda x: f"{x['Position']}:{x['Ref']}>{x['Alt']}", axis=1)
        
        # Count shared mutations
        mutation_country = data.groupby(["Mutation", "Type"])["Country"].apply(list).reset_index()
        mutation_country["Country_Count"] = mutation_country["Country"].apply(len)
        mutation_country["Countries"] = mutation_country["Country"].apply(lambda x: ", ".join(sorted(set(x))))
        
        # Shared mutations
        st.subheader("Shared Mutations Across Countries")
        shared = mutation_country[mutation_country["Country_Count"] > 1].sort_values("Country_Count", ascending=False)
        
        if not shared.empty:
            fig = px.bar(shared.head(20), x="Mutation", y="Country_Count", color="Type",
                        hover_data=["Countries"], labels={"Country_Count": "Number of Countries"})
            st.plotly_chart(fig, use_container_width=True)
        else:
            st.info("No mutations are shared across multiple countries.")
        
        # Country similarity matrix based on shared mutations
        st.subheader("Country Similarity Matrix")
        
        countries = sorted(data["Country"].unique())
        similarity_matrix = np.zeros((len(countries), len(countries)))
        
        # Build similarity matrix
        country_mutations = {country: set(data[data["Country"] == country]["Mutation"]) 
                            for country in countries}
        
        for i, country1 in enumerate(countries):
            for j, country2 in enumerate(countries):
                if i == j:
                    similarity_matrix[i, j] = 1.0  # Self-similarity
                else:
                    # Jaccard similarity
                    intersection = len(country_mutations[country1].intersection(country_mutations[country2]))
                    union = len(country_mutations[country1].union(country_mutations[country2]))
                    similarity_matrix[i, j] = intersection / union if union > 0 else 0
        
        fig = px.imshow(similarity_matrix, x=countries, y=countries, color_continuous_scale="Viridis",
                      labels=dict(x="Country", y="Country", color="Similarity"))
        st.plotly_chart(fig, use_container_width=True)

# Instructions on how to run the dashboard
st.sidebar.markdown("""
## How to Run

1. First run the analysis pipeline:
   ```
   python run_pipeline.py
   ```

2. Then launch this dashboard:
   ```
   pip install streamlit
   streamlit run dashboard.py
   ```
""")