#include <vulkan/vulkan.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define BINDING_COUNT 19
#define ERR_CAP 1024

typedef struct {
    char device_name[VK_MAX_PHYSICAL_DEVICE_NAME_SIZE];
    uint32_t vendor_id, device_id, device_type, api_version, driver_version;
    uint32_t queue_family, supports_float64;
    uint32_t max_invocations, max_size_x;
} PbufDeviceInfo;

typedef struct {
    VkInstance instance;
    VkPhysicalDevice physical;
    VkDevice device;
    VkQueue queue;
    uint32_t queue_family, workgroup_size;
    VkDescriptorSetLayout descriptor_layout;
    VkPipelineLayout pipeline_layout;
    VkPipeline pipeline;
    VkShaderModule shader;
} Runtime;

typedef struct { VkBuffer buffer; VkDeviceMemory memory; VkDeviceSize size; } Buffer;
typedef struct { VkPhysicalDevice dev; PbufDeviceInfo info; } Candidate;

static void fail(char *err, const char *msg, VkResult result) {
    if (err) snprintf(err, ERR_CAP, "%s (VkResult=%d)", msg, (int)result);
}

static VkResult make_instance(VkInstance *out) {
    VkApplicationInfo app = { .sType = VK_STRUCTURE_TYPE_APPLICATION_INFO,
        .pApplicationName = "pbuf-wl-vulkan", .applicationVersion = 1,
        .pEngineName = "pbuf", .engineVersion = 1, .apiVersion = VK_API_VERSION_1_1 };
    VkInstanceCreateInfo ci = { .sType = VK_STRUCTURE_TYPE_INSTANCE_CREATE_INFO,
        .pApplicationInfo = &app };
    return vkCreateInstance(&ci, NULL, out);
}

static int compare_candidate(const void *a, const void *b) {
    const Candidate *x = a, *y = b;
    int xd = x->info.device_type == VK_PHYSICAL_DEVICE_TYPE_DISCRETE_GPU;
    int yd = y->info.device_type == VK_PHYSICAL_DEVICE_TYPE_DISCRETE_GPU;
    if (xd != yd) return yd - xd;
    if (x->info.vendor_id != y->info.vendor_id) return x->info.vendor_id < y->info.vendor_id ? -1 : 1;
    if (x->info.device_id != y->info.device_id) return x->info.device_id < y->info.device_id ? -1 : 1;
    return strcmp(x->info.device_name, y->info.device_name);
}

static int candidates(VkInstance instance, Candidate **out, uint32_t *count, char *err) {
    uint32_t n = 0; VkResult r = vkEnumeratePhysicalDevices(instance, &n, NULL);
    if (r != VK_SUCCESS || n == 0) { fail(err, "no Vulkan physical devices", r); return -1; }
    VkPhysicalDevice *devices = calloc(n, sizeof(*devices));
    Candidate *valid = calloc(n, sizeof(*valid));
    if (!devices || !valid) { snprintf(err, ERR_CAP, "host allocation failed"); free(devices); free(valid); return -1; }
    r = vkEnumeratePhysicalDevices(instance, &n, devices);
    if (r != VK_SUCCESS) { fail(err, "physical-device enumeration failed", r); free(devices); free(valid); return -1; }
    uint32_t used = 0;
    for (uint32_t i=0; i<n; ++i) {
        VkPhysicalDeviceProperties p; VkPhysicalDeviceFeatures f;
        vkGetPhysicalDeviceProperties(devices[i], &p); vkGetPhysicalDeviceFeatures(devices[i], &f);
        uint32_t qn=0; vkGetPhysicalDeviceQueueFamilyProperties(devices[i], &qn, NULL);
        VkQueueFamilyProperties *qp = calloc(qn, sizeof(*qp));
        vkGetPhysicalDeviceQueueFamilyProperties(devices[i], &qn, qp);
        uint32_t qi = UINT32_MAX;
        for (uint32_t j=0; j<qn; ++j) if (qp[j].queueCount && (qp[j].queueFlags & VK_QUEUE_COMPUTE_BIT)) { qi=j; break; }
        free(qp);
        if (qi == UINT32_MAX || !f.shaderFloat64) continue;
        Candidate *c = &valid[used++]; c->dev=devices[i];
        strncpy(c->info.device_name, p.deviceName, sizeof(c->info.device_name)-1);
        c->info.vendor_id=p.vendorID; c->info.device_id=p.deviceID; c->info.device_type=p.deviceType;
        c->info.api_version=p.apiVersion; c->info.driver_version=p.driverVersion;
        c->info.queue_family=qi; c->info.supports_float64=f.shaderFloat64;
        c->info.max_invocations=p.limits.maxComputeWorkGroupInvocations;
        c->info.max_size_x=p.limits.maxComputeWorkGroupSize[0];
    }
    free(devices);
    if (!used) { snprintf(err, ERR_CAP, "no compute-capable Vulkan device with shader float64 support"); free(valid); return -1; }
    qsort(valid, used, sizeof(*valid), compare_candidate);
    *out=valid; *count=used; return 0;
}

int pbuf_vk_discover(int requested_index, PbufDeviceInfo *info, char *err) {
    VkInstance instance=VK_NULL_HANDLE; VkResult r=make_instance(&instance);
    if (r != VK_SUCCESS) { fail(err, "Vulkan instance creation failed", r); return -1; }
    Candidate *list=NULL; uint32_t n=0; int rc=candidates(instance, &list, &n, err);
    if (!rc) {
        if (requested_index < 0) requested_index=0;
        if ((uint32_t)requested_index >= n) { snprintf(err, ERR_CAP, "PBUF_VULKAN_DEVICE_INDEX=%d is invalid; valid device count is %u", requested_index, n); rc=-1; }
        else *info=list[requested_index].info;
    }
    free(list); vkDestroyInstance(instance, NULL); return rc;
}

static uint32_t memory_type(Runtime *rt, uint32_t bits, VkMemoryPropertyFlags wanted) {
    VkPhysicalDeviceMemoryProperties m; vkGetPhysicalDeviceMemoryProperties(rt->physical, &m);
    for (uint32_t i=0; i<m.memoryTypeCount; ++i)
        if ((bits & (1u<<i)) && (m.memoryTypes[i].propertyFlags & wanted) == wanted) return i;
    return UINT32_MAX;
}

static int make_buffer(Runtime *rt, VkDeviceSize size, Buffer *b, char *err) {
    b->size=size;
    VkBufferCreateInfo ci={.sType=VK_STRUCTURE_TYPE_BUFFER_CREATE_INFO,.size=size,
        .usage=VK_BUFFER_USAGE_STORAGE_BUFFER_BIT,.sharingMode=VK_SHARING_MODE_EXCLUSIVE};
    VkResult r=vkCreateBuffer(rt->device,&ci,NULL,&b->buffer);
    if(r!=VK_SUCCESS){fail(err,"buffer creation failed",r);return -1;}
    VkMemoryRequirements req; vkGetBufferMemoryRequirements(rt->device,b->buffer,&req);
    uint32_t mt=memory_type(rt,req.memoryTypeBits,VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT|VK_MEMORY_PROPERTY_HOST_COHERENT_BIT);
    if(mt==UINT32_MAX){snprintf(err,ERR_CAP,"no host-visible coherent Vulkan memory type");return -1;}
    VkMemoryAllocateInfo ai={.sType=VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO,.allocationSize=req.size,.memoryTypeIndex=mt};
    r=vkAllocateMemory(rt->device,&ai,NULL,&b->memory); if(r!=VK_SUCCESS){fail(err,"buffer memory allocation failed",r);return -1;}
    r=vkBindBufferMemory(rt->device,b->buffer,b->memory,0); if(r!=VK_SUCCESS){fail(err,"buffer bind failed",r);return -1;}
    return 0;
}

static int copy_in(Runtime *rt, Buffer *b, const void *src, char *err) {
    void *p=NULL; VkResult r=vkMapMemory(rt->device,b->memory,0,b->size,0,&p);
    if(r!=VK_SUCCESS){fail(err,"buffer map failed",r);return -1;} memcpy(p,src,(size_t)b->size); vkUnmapMemory(rt->device,b->memory); return 0;
}
static unsigned char *read_file(const char *path, size_t *size, char *err) {
    FILE *f=fopen(path,"rb"); if(!f){snprintf(err,ERR_CAP,"cannot open SPIR-V file: %s",path);return NULL;}
    fseek(f,0,SEEK_END); long n=ftell(f); rewind(f); unsigned char *p=malloc((size_t)n);
    if(!p || fread(p,1,(size_t)n,f)!=(size_t)n){snprintf(err,ERR_CAP,"cannot read SPIR-V file: %s",path);free(p);fclose(f);return NULL;}
    fclose(f);*size=(size_t)n;return p;
}

void *pbuf_vk_create(const char *spv_path, int requested_index, uint32_t requested_workgroup, PbufDeviceInfo *info, char *err) {
    Runtime *rt=calloc(1,sizeof(*rt)); if(!rt){snprintf(err,ERR_CAP,"runtime allocation failed");return NULL;}
    VkResult r=make_instance(&rt->instance); if(r!=VK_SUCCESS){fail(err,"Vulkan instance creation failed",r);free(rt);return NULL;}
    Candidate *list=NULL;uint32_t n=0;if(candidates(rt->instance,&list,&n,err))goto bad;
    if(requested_index<0)requested_index=0;
    if((uint32_t)requested_index>=n){snprintf(err,ERR_CAP,"PBUF_VULKAN_DEVICE_INDEX=%d is invalid; valid device count is %u",requested_index,n);free(list);goto bad;}
    rt->physical=list[requested_index].dev;rt->queue_family=list[requested_index].info.queue_family;*info=list[requested_index].info;free(list);
    uint32_t limit=info->max_invocations<info->max_size_x?info->max_invocations:info->max_size_x;
    uint32_t wg=requested_workgroup?requested_workgroup:256; if(wg>256)wg=256; while(wg>limit)wg>>=1; if(!wg){snprintf(err,ERR_CAP,"device has no usable compute workgroup size");goto bad;} rt->workgroup_size=wg;
    float priority=1.0f;VkDeviceQueueCreateInfo qci={.sType=VK_STRUCTURE_TYPE_DEVICE_QUEUE_CREATE_INFO,.queueFamilyIndex=rt->queue_family,.queueCount=1,.pQueuePriorities=&priority};
    VkPhysicalDeviceFeatures features={.shaderFloat64=VK_TRUE};VkDeviceCreateInfo dci={.sType=VK_STRUCTURE_TYPE_DEVICE_CREATE_INFO,.queueCreateInfoCount=1,.pQueueCreateInfos=&qci,.pEnabledFeatures=&features};
    r=vkCreateDevice(rt->physical,&dci,NULL,&rt->device);if(r!=VK_SUCCESS){fail(err,"logical device creation failed",r);goto bad;}vkGetDeviceQueue(rt->device,rt->queue_family,0,&rt->queue);
    VkDescriptorSetLayoutBinding bindings[BINDING_COUNT];memset(bindings,0,sizeof(bindings));
    for(uint32_t i=0;i<BINDING_COUNT;++i){bindings[i].binding=i;bindings[i].descriptorType=VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;bindings[i].descriptorCount=1;bindings[i].stageFlags=VK_SHADER_STAGE_COMPUTE_BIT;}
    VkDescriptorSetLayoutCreateInfo dl={.sType=VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO,.bindingCount=BINDING_COUNT,.pBindings=bindings};
    r=vkCreateDescriptorSetLayout(rt->device,&dl,NULL,&rt->descriptor_layout);if(r!=VK_SUCCESS){fail(err,"descriptor layout creation failed",r);goto bad;}
    VkPipelineLayoutCreateInfo pl={.sType=VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO,.setLayoutCount=1,.pSetLayouts=&rt->descriptor_layout};
    r=vkCreatePipelineLayout(rt->device,&pl,NULL,&rt->pipeline_layout);if(r!=VK_SUCCESS){fail(err,"pipeline layout creation failed",r);goto bad;}
    size_t spv_size=0;unsigned char *spv=read_file(spv_path,&spv_size,err);if(!spv)goto bad;
    VkShaderModuleCreateInfo sm={.sType=VK_STRUCTURE_TYPE_SHADER_MODULE_CREATE_INFO,.codeSize=spv_size,.pCode=(uint32_t*)spv};r=vkCreateShaderModule(rt->device,&sm,NULL,&rt->shader);free(spv);if(r!=VK_SUCCESS){fail(err,"shader module creation failed",r);goto bad;}
    VkSpecializationMapEntry map={.constantID=0,.offset=0,.size=sizeof(uint32_t)};VkSpecializationInfo spec={.mapEntryCount=1,.pMapEntries=&map,.dataSize=sizeof(uint32_t),.pData=&rt->workgroup_size};
    VkPipelineShaderStageCreateInfo stage={.sType=VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO,.stage=VK_SHADER_STAGE_COMPUTE_BIT,.module=rt->shader,.pName="main",.pSpecializationInfo=&spec};
    VkComputePipelineCreateInfo pc={.sType=VK_STRUCTURE_TYPE_COMPUTE_PIPELINE_CREATE_INFO,.stage=stage,.layout=rt->pipeline_layout};r=vkCreateComputePipelines(rt->device,VK_NULL_HANDLE,1,&pc,NULL,&rt->pipeline);if(r!=VK_SUCCESS){fail(err,"compute pipeline creation failed",r);goto bad;}
    return rt;
bad:
    if(rt->device){if(rt->shader)vkDestroyShaderModule(rt->device,rt->shader,NULL);if(rt->pipeline_layout)vkDestroyPipelineLayout(rt->device,rt->pipeline_layout,NULL);if(rt->descriptor_layout)vkDestroyDescriptorSetLayout(rt->device,rt->descriptor_layout,NULL);vkDestroyDevice(rt->device,NULL);}if(rt->instance)vkDestroyInstance(rt->instance,NULL);free(rt);return NULL;
}

void pbuf_vk_destroy(void *handle) {
    Runtime *rt=handle;if(!rt)return;vkDeviceWaitIdle(rt->device);vkDestroyPipeline(rt->device,rt->pipeline,NULL);vkDestroyShaderModule(rt->device,rt->shader,NULL);vkDestroyPipelineLayout(rt->device,rt->pipeline_layout,NULL);vkDestroyDescriptorSetLayout(rt->device,rt->descriptor_layout,NULL);vkDestroyDevice(rt->device,NULL);vkDestroyInstance(rt->instance,NULL);free(rt);
}

int pbuf_vk_propagate(void *handle, const double **inputs, const uint64_t *lengths,
                      double **outputs, uint32_t ray_count, uint32_t nx, uint32_t ny,
                      double step, uint32_t steps, const double *checkpoints,
                      uint32_t checkpoint_count, char *err) {
    Runtime *rt=handle; Buffer b[BINDING_COUNT];memset(b,0,sizeof(b));
    VkDescriptorPool pool=VK_NULL_HANDLE;VkDescriptorSet set=VK_NULL_HANDLE;VkCommandPool cp=VK_NULL_HANDLE;VkCommandBuffer cmd=VK_NULL_HANDLE;VkFence fence=VK_NULL_HANDLE;int rc=-1;
    double config[6]={(double)ray_count,(double)nx,(double)ny,step,(double)steps,(double)checkpoint_count};
    VkDeviceSize sizes[BINDING_COUNT];
    for(int i=0;i<10;++i)sizes[i]=lengths[i]*sizeof(double);
    for(int i=10;i<16;++i)sizes[i]=(VkDeviceSize)ray_count*checkpoint_count*sizeof(double);
    sizes[16]=sizeof(config);sizes[17]=(VkDeviceSize)checkpoint_count*sizeof(double);
    sizes[18]=(VkDeviceSize)ray_count*sizeof(double);
    for(int i=0;i<BINDING_COUNT;++i)if(make_buffer(rt,sizes[i],&b[i],err))goto done;
    for(int i=0;i<10;++i)if(copy_in(rt,&b[i],inputs[i],err))goto done;
    if(copy_in(rt,&b[16],config,err)||copy_in(rt,&b[17],checkpoints,err))goto done;
    VkDescriptorPoolSize ps={.type=VK_DESCRIPTOR_TYPE_STORAGE_BUFFER,.descriptorCount=BINDING_COUNT};VkDescriptorPoolCreateInfo pci={.sType=VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO,.maxSets=1,.poolSizeCount=1,.pPoolSizes=&ps};VkResult r=vkCreateDescriptorPool(rt->device,&pci,NULL,&pool);if(r!=VK_SUCCESS){fail(err,"descriptor pool creation failed",r);goto done;}
    VkDescriptorSetAllocateInfo sai={.sType=VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO,.descriptorPool=pool,.descriptorSetCount=1,.pSetLayouts=&rt->descriptor_layout};r=vkAllocateDescriptorSets(rt->device,&sai,&set);if(r!=VK_SUCCESS){fail(err,"descriptor set allocation failed",r);goto done;}
    VkDescriptorBufferInfo bi[BINDING_COUNT];VkWriteDescriptorSet wr[BINDING_COUNT];memset(wr,0,sizeof(wr));
    for(uint32_t i=0;i<BINDING_COUNT;++i){bi[i]=(VkDescriptorBufferInfo){.buffer=b[i].buffer,.offset=0,.range=b[i].size};wr[i]=(VkWriteDescriptorSet){.sType=VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET,.dstSet=set,.dstBinding=i,.descriptorCount=1,.descriptorType=VK_DESCRIPTOR_TYPE_STORAGE_BUFFER,.pBufferInfo=&bi[i]};}vkUpdateDescriptorSets(rt->device,BINDING_COUNT,wr,0,NULL);
    VkCommandPoolCreateInfo cpi={.sType=VK_STRUCTURE_TYPE_COMMAND_POOL_CREATE_INFO,.queueFamilyIndex=rt->queue_family};r=vkCreateCommandPool(rt->device,&cpi,NULL,&cp);if(r!=VK_SUCCESS){fail(err,"command pool creation failed",r);goto done;}
    VkCommandBufferAllocateInfo cai={.sType=VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO,.commandPool=cp,.level=VK_COMMAND_BUFFER_LEVEL_PRIMARY,.commandBufferCount=1};r=vkAllocateCommandBuffers(rt->device,&cai,&cmd);if(r!=VK_SUCCESS){fail(err,"command buffer allocation failed",r);goto done;}
    VkCommandBufferBeginInfo cb={.sType=VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO,.flags=VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT};vkBeginCommandBuffer(cmd,&cb);vkCmdBindPipeline(cmd,VK_PIPELINE_BIND_POINT_COMPUTE,rt->pipeline);vkCmdBindDescriptorSets(cmd,VK_PIPELINE_BIND_POINT_COMPUTE,rt->pipeline_layout,0,1,&set,0,NULL);vkCmdDispatch(cmd,(ray_count+rt->workgroup_size-1)/rt->workgroup_size,1,1);r=vkEndCommandBuffer(cmd);if(r!=VK_SUCCESS){fail(err,"command buffer recording failed",r);goto done;}
    VkFenceCreateInfo fi={.sType=VK_STRUCTURE_TYPE_FENCE_CREATE_INFO};r=vkCreateFence(rt->device,&fi,NULL,&fence);if(r!=VK_SUCCESS){fail(err,"fence creation failed",r);goto done;}VkSubmitInfo si={.sType=VK_STRUCTURE_TYPE_SUBMIT_INFO,.commandBufferCount=1,.pCommandBuffers=&cmd};r=vkQueueSubmit(rt->queue,1,&si,fence);if(r==VK_SUCCESS)r=vkWaitForFences(rt->device,1,&fence,VK_TRUE,UINT64_MAX);if(r!=VK_SUCCESS){fail(err,"compute dispatch failed",r);goto done;}
    for(int i=10;i<16;++i){void *p=NULL;r=vkMapMemory(rt->device,b[i].memory,0,b[i].size,0,&p);if(r!=VK_SUCCESS){fail(err,"output map failed",r);goto done;}memcpy(outputs[i-10],p,(size_t)b[i].size);vkUnmapMemory(rt->device,b[i].memory);}
    { void *p=NULL;r=vkMapMemory(rt->device,b[18].memory,0,b[18].size,0,&p);if(r!=VK_SUCCESS){fail(err,"unit-error output map failed",r);goto done;}memcpy(outputs[6],p,(size_t)b[18].size);vkUnmapMemory(rt->device,b[18].memory); }
    rc=0;
done:
    if (fence) vkDestroyFence(rt->device, fence, NULL);
    if (cp) vkDestroyCommandPool(rt->device, cp, NULL);
    if (pool) vkDestroyDescriptorPool(rt->device, pool, NULL);
    for (int i=0; i<BINDING_COUNT; ++i) {
        if (b[i].buffer) vkDestroyBuffer(rt->device, b[i].buffer, NULL);
        if (b[i].memory) vkFreeMemory(rt->device, b[i].memory, NULL);
    }
    return rc;
}

uint32_t pbuf_vk_workgroup_size(void *handle) { return ((Runtime*)handle)->workgroup_size; }

int pbuf_vk_kde(void *handle, const double *u, const double *v, double *output,
                uint32_t n, const double *bandwidth, char *err) {
    Runtime *rt=handle; Buffer b[BINDING_COUNT];memset(b,0,sizeof(b));
    VkDescriptorPool pool=VK_NULL_HANDLE;VkCommandPool cp=VK_NULL_HANDLE;VkCommandBuffer cmd=VK_NULL_HANDLE;VkFence fence=VK_NULL_HANDLE;int rc=-1;
    double config[5]={(double)n,bandwidth[0],bandwidth[1],0.0,0.0}; VkDeviceSize sizes[BINDING_COUNT];
    for(int i=0;i<BINDING_COUNT;++i)sizes[i]=sizeof(double);
    sizes[0]=sizes[1]=sizes[2]=(VkDeviceSize)n*sizeof(double);sizes[3]=sizeof(config);
    for(int i=0;i<BINDING_COUNT;++i)if(make_buffer(rt,sizes[i],&b[i],err))goto done;
    if(copy_in(rt,&b[0],u,err)||copy_in(rt,&b[1],v,err)||copy_in(rt,&b[3],config,err))goto done;
    VkDescriptorPoolSize ps={.type=VK_DESCRIPTOR_TYPE_STORAGE_BUFFER,.descriptorCount=BINDING_COUNT};
    VkDescriptorPoolCreateInfo pci={.sType=VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO,.maxSets=1,.poolSizeCount=1,.pPoolSizes=&ps};VkResult r=vkCreateDescriptorPool(rt->device,&pci,NULL,&pool);if(r!=VK_SUCCESS){fail(err,"KDE descriptor pool failed",r);goto done;}
    VkDescriptorSet set=VK_NULL_HANDLE;VkDescriptorSetAllocateInfo sai={.sType=VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO,.descriptorPool=pool,.descriptorSetCount=1,.pSetLayouts=&rt->descriptor_layout};r=vkAllocateDescriptorSets(rt->device,&sai,&set);if(r!=VK_SUCCESS){fail(err,"KDE descriptor set failed",r);goto done;}
    VkDescriptorBufferInfo bi[BINDING_COUNT];VkWriteDescriptorSet wr[BINDING_COUNT];memset(wr,0,sizeof(wr));
    for(uint32_t i=0;i<BINDING_COUNT;++i){bi[i]=(VkDescriptorBufferInfo){.buffer=b[i].buffer,.offset=0,.range=b[i].size};wr[i]=(VkWriteDescriptorSet){.sType=VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET,.dstSet=set,.dstBinding=i,.descriptorCount=1,.descriptorType=VK_DESCRIPTOR_TYPE_STORAGE_BUFFER,.pBufferInfo=&bi[i]};}vkUpdateDescriptorSets(rt->device,BINDING_COUNT,wr,0,NULL);
    VkCommandPoolCreateInfo cpi={.sType=VK_STRUCTURE_TYPE_COMMAND_POOL_CREATE_INFO,.queueFamilyIndex=rt->queue_family};r=vkCreateCommandPool(rt->device,&cpi,NULL,&cp);if(r!=VK_SUCCESS){fail(err,"KDE command pool failed",r);goto done;}
    VkCommandBufferAllocateInfo cai={.sType=VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO,.commandPool=cp,.level=VK_COMMAND_BUFFER_LEVEL_PRIMARY,.commandBufferCount=1};vkAllocateCommandBuffers(rt->device,&cai,&cmd);
    VkFenceCreateInfo fi={.sType=VK_STRUCTURE_TYPE_FENCE_CREATE_INFO};vkCreateFence(rt->device,&fi,NULL,&fence);
    for(uint32_t offset=0;offset<n;offset+=1024u){
      uint32_t count=n-offset<1024u?n-offset:1024u;config[3]=(double)offset;config[4]=(double)count;
      if(copy_in(rt,&b[3],config,err))goto done;
      vkResetCommandPool(rt->device,cp,0);VkCommandBufferBeginInfo cb={.sType=VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO,.flags=VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT};vkBeginCommandBuffer(cmd,&cb);vkCmdBindPipeline(cmd,VK_PIPELINE_BIND_POINT_COMPUTE,rt->pipeline);vkCmdBindDescriptorSets(cmd,VK_PIPELINE_BIND_POINT_COMPUTE,rt->pipeline_layout,0,1,&set,0,NULL);vkCmdDispatch(cmd,(count+rt->workgroup_size-1)/rt->workgroup_size,1,1);vkEndCommandBuffer(cmd);
      vkResetFences(rt->device,1,&fence);VkSubmitInfo si={.sType=VK_STRUCTURE_TYPE_SUBMIT_INFO,.commandBufferCount=1,.pCommandBuffers=&cmd};r=vkQueueSubmit(rt->queue,1,&si,fence);if(r==VK_SUCCESS)r=vkWaitForFences(rt->device,1,&fence,VK_TRUE,UINT64_MAX);if(r!=VK_SUCCESS){fail(err,"KDE dispatch failed",r);goto done;}
    }
    {void *p=NULL;r=vkMapMemory(rt->device,b[2].memory,0,b[2].size,0,&p);if(r!=VK_SUCCESS){fail(err,"KDE output map failed",r);goto done;}memcpy(output,p,(size_t)b[2].size);vkUnmapMemory(rt->device,b[2].memory);}rc=0;
done:
    if(fence)vkDestroyFence(rt->device,fence,NULL);if(cp)vkDestroyCommandPool(rt->device,cp,NULL);if(pool)vkDestroyDescriptorPool(rt->device,pool,NULL);
    for(int i=0;i<BINDING_COUNT;++i){if(b[i].buffer)vkDestroyBuffer(rt->device,b[i].buffer,NULL);if(b[i].memory)vkFreeMemory(rt->device,b[i].memory,NULL);}return rc;
}
