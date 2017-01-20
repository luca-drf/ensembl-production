-- Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
-- Copyright [2016-2017] EMBL-European Bioinformatics Institute
-- 
-- Licensed under the Apache License, Version 2.0 (the "License");
-- you may not use this file except in compliance with the License.
-- You may obtain a copy of the License at
-- 
--      http://www.apache.org/licenses/LICENSE-2.0
-- 
-- Unless required by applicable law or agreed to in writing, software
-- distributed under the License is distributed on an "AS IS" BASIS,
-- WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
-- See the License for the specific language governing permissions and
-- limitations under the License.

ALTER TABLE biotype MODIFY COLUMN db_type set('cdna', 'core', 'coreexpressionatlas','coreexpressionest', 'coreexpressiongnf', 'funcgen',
                                              'otherfeatures', 'rnaseq', 'variation', 'vega', 'presite', 'sangerverga') NOT NULL DEFAULT 'core';
ALTER TABLE biotype ADD COLUMN biotype_group ENUM('coding','pseudogene','snoncoding','lnoncoding','undefined') DEFAULT NULL;

-- Patch identifier
INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'patch', 'patch_72_73c.sql|biotype group in biotype table');
