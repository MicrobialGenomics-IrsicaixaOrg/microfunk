# Function fetch_humann_db() @ L155

    Code
      fetch_humann_db("test")
    Message
      
      -- Fetching HUMAnN3 databases --------------------------------------------------
      i Cache directory .microfunk already exists
      v File test_v201901b.tsv.gz cached

---

    Code
      fetch_humann_db("test")
    Message
      
      -- Fetching HUMAnN3 databases --------------------------------------------------
      i Cache directory .microfunk already exists
      v File test_v201901b.tsv.gz cached

---

    Code
      fetch_humann_db("test", overwrite = TRUE)
    Message
      
      -- Fetching HUMAnN3 databases --------------------------------------------------
      i Cache directory .microfunk already exists
      v File test_v201901b.tsv.gz cached
      i Fetching test_v201901b.tsv.gz from AWS
      v File test_v201901b.tsv.gz downloaded to ~/.microfunk

